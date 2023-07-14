#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Alexis Perez-Bellido
# Code for wm task. This is use to test healthy participants. 
# There are some catch events to maintain eyes on fixation point during trial (only used in healthy participants)


#pip install -U pip setuptools
#pip install tobii_research

version = 0.1

import numpy as np
from scipy.stats import vonmises
import matplotlib.pyplot as plt
import os, sys, ctypes, csv
import exp_func as exp
import stimuli as st
import instructions as instr
import serial
from time import sleep,time
import datetime
from psychopy.hardware import keyboard



from random import random
from numpy import sin, pi, cos
from scipy import signal, stats
from numpy.linalg import norm
from pandas import DataFrame


from psychopy import visual, logging, core, event,  gui, data
from psychopy.tools.filetools import fromFile, toFile # wrappers to save pickles
from psychopy.preferences import prefs
from psychopy import prefs

import importlib
#importlib.reload(exp)

#prefs.general['winType'] =  'glfw'
#print(prefs)

Clock = core.Clock()

# Some general presets
event.globalKeys.clear() # implementing a global event to quit the program at any time by pressing ctrl+q
event.globalKeys.add(key='q', modifiers=['ctrl'], func=core.quit)



log_on=True
sst = False # Using parallel port to send triggers
tobii = False
showEyes = False

# Initiating eyetracker
if tobii:
    import tobii_research as tr
    #import time

    found_eyetrackers = tr.find_all_eyetrackers()
    my_eyetracker = found_eyetrackers[0]
    print("Address: " + my_eyetracker.address)
    print("Model: " + my_eyetracker.model)
    print("Name (It's OK if this is empty): " + my_eyetracker.device_name)
    print("Serial number: " + my_eyetracker.serial_number)

''' Variable names
columns = ['device_time_stamp', 'system_time_stamp', 'left_gaze_direction_unit_vector', 
           'left_gaze_direction_validity', 'left_gaze_origin_position_in_hmd_coordinates', 
           'left_gaze_origin_validity', 'left_pupil_diameter', 'left_pupil_validity', 
           'left_pupil_position_in_tracking_area', 'left_pupil_position_validity', 
           'right_gaze_direction_unit_vector', 'right_gaze_direction_validity', 
           'right_gaze_origin_position_in_hmd_coordinates', 'right_gaze_origin_validity',
           'right_pupil_diameter', 'right_pupil_validity', 'right_pupil_position_in_tracking_area', 
           'right_pupil_position_validity']
'''

global_gaze_data = {}
global_gaze_data['left_pos'] = np.array([np.nan,np.nan])
global_gaze_data['right_pos'] = np.array([np.nan,np.nan])
global_gaze_data['left_pupil'] = np.nan
global_gaze_data['right_pupil'] = np.nan
global_gaze_data['left_pupil_validity'] = np.nan
global_gaze_data['right_pupil_validity'] = np.nan
global_gaze_data['device_time'] = np.nan
global_gaze_data['system_time'] = np.nan

all_gaze_data = [] # saving all the eye data here


if tobii:
    
    def gaze_data_callback(gaze_data): # this is a callback function that will return values
        global global_gaze_data
        global_gaze_data['left_pos'] = gaze_data['left_gaze_point_on_display_area']
        global_gaze_data['right_pos'] = gaze_data['right_gaze_point_on_display_area']
        global_gaze_data['left_pupil'] = gaze_data['left_pupil_diameter']
        global_gaze_data['right_pupil'] = gaze_data['right_pupil_diameter']
        global_gaze_data['left_pupil_validity'] = gaze_data['left_pupil_validity']
        global_gaze_data['right_pupil_validity'] = gaze_data['right_pupil_validity']
        global_gaze_data['device_time'] = gaze_data['device_time_stamp']
        global_gaze_data['system_time'] = gaze_data['system_time_stamp']
        all_gaze_data.append(gaze_data)



    # function to combine eye and left eye data
    def get_combined_eyes(gdata):
        combined_eyes = {}
        LPos = np.array(gdata['left_pos'])
        RPos = np.array(gdata['right_pos'])
        combined_eyes['EyesPos'] = np.nanmean([LPos,RPos], axis = 0)

        LPup = np.array(gdata['left_pupil'])
        RPup = np.array(gdata['right_pupil'])
        combined_eyes['EyesPup'] = np.nanmean([LPup,RPup], axis = 0)  
        return combined_eyes

#Triggers

tg_fp =  '07'
tg_stim = '27'
tg_block = '61'
tg_resphase = '41'
tg_respinit = '35'
tg_resp = '31'
tg_zero = '00'

# Collect subject data
expInfo, resultspath = exp.mainexp_subject_info(version)
subj_id = expInfo['subjInfo']['observer']

# Loading monitor definitions
monitores = st.monitor_def() 

mon, expInfo['monitor'] = exp.define_monitor(monitores[4]) # select the correct monitor

# Creating a new experimental window
monitor_features = {}
monitor_features['monitor'] = mon
monitor_features['units'] = 'cm' # units to define your stimuli
monitor_features['screen_id'] = 0 # 1 when using a extended display 
monitor_features['full']  = True
monitor_features['Hz'] =  60 #'auto' #144 #60 this can be set to "auto" to estimate the refreshing rate of the monitor, although it can fail often
res =  expInfo['monitor']['monitor_pixels']
monitor_features['resolution'] = res
kb = keyboard.Keyboard()
   
win, monitor_features = exp.create_window(monitor_features)
ifi = monitor_features['ifi']


# Parameters
radius = 200  # Radius of the circular grating
sf = 0.05  # Spatial frequency of the grating
ori = 45  # Orientation of the grating
mask_radius = 250  # Radius of the circular mask

st_prop = st.stim_config(ifi) #Loading stim properties

fixation = visual.PatchStim(win=win, mask='gauss', size=1, pos=[0,0], sf=0, color=[-1,-1,-1],units='cm')
fixation_ct = visual.PatchStim(win=win, mask='gauss', size=1, pos=[0,0], sf=0, color=[-0.75,-0.75,-0.75],units='cm')


resp = visual.Circle(win,radius=0.7,pos=[0,0],fillColor='red',units='cm')
ring = visual.Circle(win,radius=st_prop['radius'],pos=[0,0],edges = 60, lineColor='black',units='cm', fillColor = [0,0,0], interpolate = True)


stim_HC = visual.Circle(win,radius=st_prop['size_CUE'] ,pos=[0,0],fillColor=[1,1,1],units= "deg")
stim_LC = visual.Circle(win,radius=st_prop['size_CUE'],pos=[0,0],fillColor=[0.45,0.45,0.45],units= "deg")

# Circular grating
radialGrating = visual.RadialStim(win, tex='sinXsin', color=0.45, size=st_prop['size_CUE']*2, radialCycles=4, angularCycles=0,
    autoLog=False,units= "deg")  # this stim changes too much for autologging to be useful



# Creating the stimuli positions
stims_pos_angle = np.linspace(0,360,11)[0:-1] # spacing the circular space in 10 positions
stims_pos_angle = stims_pos_angle + np.random.uniform(-22.5,22.5,1) # adding some noise to the positions
stims_pos_angle = np.random.permutation(stims_pos_angle) # randomize the order of the positions

# Creating the trial sequences of orientations for this participant based on randomly generated pairs
# random pair combinations
npairs = int(np.shape(stims_pos_angle)[0] / 2) 
pairs_combs = np.reshape(stims_pos_angle, (npairs, 2)) # creating pairs of positions in 5 combinations
#pairs_combs = np.tile(pairs_combs, (22,1))
nTrials = 20 # number of trials
nblocks = 8

if subj_id == 'test':
    nTrials = 10 # number of trials
    nblocks = 1


mouse=event.Mouse(win=win,visible=False)

# output files paths (csv and psychopy pickle)
date = datetime.datetime.now().strftime('%d%m%Y%H%M')
if subj_id == 'test':
    filename = subj_id + '.csv'
    filepy = subj_id + '.psydat'
    filepeye = subj_id + '_eye.psydat'
else:
    filename = subj_id + '_' + date + '.csv'
    filepy = subj_id + '.psydat'
    filepeye = subj_id + '_eye.psydat'

filename = resultspath +'/'+ filename
filepy = resultspath +'/'+ filepy
filepeye = resultspath +'/'+ filepeye



# variables recorded
cnames = ['subj','trial','block', 'delay','fix','trial_seq',
'last_stim_pairpos','dev_ID', 'dev_angle', 'dev_pairpos' , 'dev_times',  'devRT','keypressed',
'radious','T_Angle','resp_mod',
'choice_x','choice_y',  'choiceAngle','choiceR','wmRT','movTime',
'm_pos_x','m_pos_y','m_clock']



if log_on:
    with open(filename, 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=';',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerow(cnames)
            


ct_frames = np.round(50/ifi) # duration of catch trial

def routines(dim):
    exp.quit_handler(win)      
    if dim: 
        fixation_ct.draw()
    else:
            fixation.draw()
    win.update()
    x_m,y_m=mouse.getPos()  


mouse = event.Mouse(visible = True)
mouse.setPos([0,0])


## Experiment design

def circdist(angles1,angles2):       
    return np.angle(np.exp(1j*angles1)/np.exp(1j*angles2))


cond_resps = ['drag']
prior_kappa = 2 # kappa for the von mises distribution

 
main_exp  = {}
main_exp['nblocks']     = nblocks # 4 # totaltime = 90 * 6 * 5
main_exp['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing block data
main_exp['Pair_Combinations'] = pairs_combs


if tobii:
    main_eye = {}
    main_eye['Exp_blocks']  = [None] * main_exp['nblocks'] # assigning memory for storing EYE block data
    xsize =   expInfo['monitor']['monitor_width']
    ysize =   expInfo['monitor']['monitor_height']


main_exp['trial_reps']  = 25 # 25 trials * 2 delays = 50 
if subj_id == 'test':
    main_exp['trial_reps']  = 8 # 8 * 2 delays = 16
    sst = False


# Initiating parallel port
if sst: 
    p_port = serial.Serial('COM3', 115200, timeout = 0 )
    p_port.write(b'00')
    core.wait(0.2)
    p_port.write(b'RR')



win.mouseVisible = False 


if subj_id == 'test':
    instr.main_instructions(win)

if tobii:
    my_eyetracker.subscribe_to(tr.EYETRACKER_GAZE_DATA, gaze_data_callback, as_dictionary=True)
    #my_eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, gaze_data_callback)


for thisBlock in range(main_exp['nblocks']): # iterate over blocks
    #thisBlock = 0 
    resp_mod = 'drag' 
    print(thisBlock)
        
    # lets trigger the beggining of the experiment
    if sst: win.callOnFlip(p_port.write, tg_block.encode())
    win.flip()
    if sst: win.callOnFlip(p_port.write, tg_zero.encode())
    win.flip()
    
    BlockClockStart = Clock.getTime() # block experiment time
    block = {} # dummy variable to save block data
    
    trialClock = core.Clock()

    trialClocktimes = np.array([]) # saving whole times here
    correct_seq = np.array([]) # saving seq. error responses per trial sequence

    trials = [None] * nTrials
    kfix = True # is the fixation correct?
                    

    # Generating the trial sequence for this block extracting pairs from the full sequence
    for itrials in range(nTrials):
        #itrials = 0
        # initialize trial vector
        thr_trials_var = [cnames]

        npairsxtrial = np.random.randint(10, 20, 1)[0]   
        # Randomly sample npairsxtrial rows from the array
        sampled_rows = np.random.choice(pairs_combs.shape[0], size= npairsxtrial, replace=True)
        # Select the sampled rows from the original array
        seq_pairs_combs =  pairs_combs[sampled_rows, :]
        trials[itrials] = seq_pairs_combs.flatten()

    
    instr.block_start(win, thisBlock, main_exp['nblocks'], resp_mod )    

    eyedata = [None] * nTrials
    

    for n_trial, thisTrial in enumerate(trials): # iterate over trials
        #n_trial = 0
        stim_oris =  trials[n_trial]      
        last_stim = np.random.choice(['leading', 'trailing'])
        if last_stim == 'leading': 
            stim_oris = stim_oris[:-1] 

        deviants = np.random.randint(2,  len(stim_oris)-4, 2)[:] # 2 deviants per tria, cannot happen in the first presentation nor in 5 before the last one
        if 5 > np.abs(np.diff(deviants)): # if the difference between the two deviants is less than 5 other stimuli, just plot one of the two
            deviants = np.random.choice(deviants)
        
        #initialize eye data variables for this trial
        if tobii:
            trial_eYe_data = {}
            eyes_stamps = []   
            eyes_stamps.append(tr.get_system_time_stamp()) # first time stamp        

        trep = True
        
        while trep:
            # during your trial
            responses = []
            deviant_times = []
            deviant_angles = []
            deviant_pos = []
           
            trialClockStart = Clock.getTime() 
            # starting trial sequence presentation
            for istim, pos_angle in enumerate(stim_oris):

                # mapping orientation in x y coordinates
                xstim,ystim = exp.toCar(st_prop['radius'] ,pos_angle)  # 1.5 corrects for the different in radius between the cue and the ring
                # plotting deviant target or distractor
                if np.isin(istim, deviants):
                    stim_ID = radialGrating
                else:
                    stim_ID = stim_LC

                stim_ID.setPos([xstim,ystim])

                if kfix == True:
                    eyepos = []
                    for frameN in range(st_prop['CUE_frames']): 
                        if frameN == 0:   
                            if np.isin(istim, deviants): 
                                # saving data about the deviant stimulus
                                deviant_angles.append(pos_angle) 

                                if np.isin(pos_angle,  pairs_combs[:,1]): # if is in the second position is a trailing
                                    deviant_pos.append('trailing')
                                else:
                                    deviant_pos.append('leading')                              

                                deviant_times.append(Clock.getTime())
                            if sst: win.callOnFlip(p_port.write, tg_stim.encode()) # send trigger
                        if frameN == 1:     
                            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
                        stim_ID.draw()
                        if tobii:
                            Eyes = get_combined_eyes(global_gaze_data)
                            eyepos.append(Eyes['EyesPos'])
                            if frameN == st_prop['CUE_frames']-1:
                                EyesPos = np.nanmean(eyepos, axis=0)
                                if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                                if  EyesPos[0] >  eye_lim or EyesPos[0] < 1-eye_lim or EyesPos[1] >  eye_lim or EyesPos[1] < 1-eye_lim:
                                    kfix = False
                                    break
                        thisResp = exp.getResponse(win, ["space"], Clock) 
                        if thisResp != None:
                            responses.append(thisResp)
                        routines(False)
                   
                    for frameN in range(st_prop['CUE_frames']):
                        thisResp = exp.getResponse(win, ["space"], Clock) 
                        if thisResp != None:
                            responses.append(thisResp)
                        routines(False)


            # get response timings
            if np.array(np.shape(responses))[0] != 0: 
                dat = np.array(responses) 
                rt = dat[:, 0].astype(float) - trialClockStart
            deviant_times = np.array(deviant_times) - trialClockStart
            

            angle = pos_angle#vonmises.rvs(10**-10, 0, size=1)[0]
            #angle = np.rad2deg(angle)
            #if angle > 180: angle = angle -360
            #if angle < -180: angle = angle +360
         

            # memory delay
            t_delay = 3000 #thisTrial['delay'] 
            delay_frames = int(t_delay/ifi)

            # is this trial a catch trial?
            ct = False
            if subj_id == 'test':
                if np.random.randint(low = 0, high = 10) < 3: # for the demo increase the number of catch trials 
                    ct = True
            else:
                if np.random.randint(low = 0, high = 10) == 0: 
                    ct = True
            ct = False # uncomment this for healthy participants
            
            # decidign when the CT should appear
            ct_frames_interval =  [(ct_frames * 2), delay_frames-(ct_frames * 2)] # cant be smaller than ct duration, and leave some time from onset and offset of delay
            delay_ct = np.random.randint(low = ct_frames_interval[0], high= ct_frames_interval[1]) 
        
            
            trialClock.reset()
            ts_d=time()
            pos = mouse.getPos()    
            # DELAY PERIOD

            if tobii:
                eyes_stamps.append(tr.get_system_time_stamp()) # first time stamp

            if kfix == True:
                eyepos = []
                for frameN in range(delay_frames):
                    if tobii:
                        Eyes = get_combined_eyes(global_gaze_data)
                        eyepos.append(Eyes['EyesPos'])
                        if frameN == delay_frames-1: # in last frame, veriify whether the eyes are still in the fixation window
                            EyesPos = np.nanmean(eyepos, axis=0)
                            if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                            if  EyesPos[0] >  eye_lim or EyesPos[0] < 1-eye_lim or EyesPos[1] >  eye_lim or EyesPos[1] < 1-eye_lim:
                                kfix = False
                                break 
                                        
                    if ct and ((frameN >= delay_ct) and (frameN < delay_ct+ct_frames)):
                        routines(True)
                    else:
                        routines(False)
                        res

            
            #mouse.setPos((0,0))
            #pos = mouse.getPos()
            choice_trial = 0
            event.clearEvents()
            #fixation.setColor("black")
            
            m_pos_x=[]; m_pos_y=[]
            m_clock=[]     
            
            trialClock.reset()
            ts_r = time()
            
            rinit_time = trialClock.getTime()                        # response time
            kb.clock.reset()  # when you want to start the timer from
            kb.clearEvents()
            keypressed = None

            mouse.setPos((0,0))
                                           
            posresp = (0, 0) # initialising the position of the response
            if tobii:
                eyes_stamps.append(tr.get_system_time_stamp()) # response init time stamp
            
            if kfix == True:
                while (mouse.getPressed()[0]==0) and (keypressed == None ): # and frameN in range(RESPONSE_MAX)):
                    keys = kb.getKeys(keyList = ['space'], waitRelease=True, clear=True)
                    for thisKey in keys:
                        keypressed = thisKey.name
                        print(keypressed, thisKey.tDown, thisKey.rt)
                        
                        
                    if resp_mod == 'drag':
                        x,y=pos=mouse.getPos()
                        posresp = (x, y) 
                        
                        t=trialClock.getTime() 
                        m_pos_x.append(x)
                        m_pos_y.append(y)
                        m_clock.append(t)
                        resp.setPos(posresp)
                        resp.draw()
                    else:
                        if tobii:
                            Eyes = get_combined_eyes(global_gaze_data)
                            #eyepos.append(Eyes['EyesPos'])
                            #if frameN == delay_frames-1: # in last frame, veriify whether the eyes are still in the fixation window
                             #   EyesPos = np.nanmean(eyepos, axis=0)
                              #  if np.isnan(EyesPos[0]):  EyesPos = np.array([0.0,0.0])
                              
                                
                        if showEyes:
                            Eyes['EyesPos'] * res
                            if  ~np.isnan(Eyes['EyesPos'][0]):
                                posresp = (Eyes['EyesPos']-0.5) *np.array([xsize,-ysize])#* res
                                resp.setPos(posresp)
                                resp.draw() 
                                
                                m_pos_x.append(posresp[0])
                                m_pos_y.append(posresp[1])
                                m_clock.append(trialClock.getTime())
                            
                   
                    if frameN == 0:     
                        if sst: win.callOnFlip(p_port.write, tg_respinit.encode()) # send trigger
                    if frameN == 1:     
                        if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins             
                    routines(False)
                    frameN += 1
                    pass
          
            #if kfix == True:
            resp_time = trialClock.getTime() - rinit_time                  # response time                       # response time
            rep_pos=posresp
                
            angle_T = angle
            r_T = st_prop['radius']
            ts_e = time()                                       # end time                   
            mtime = trialClock.getTime()                        # movement time 
            
            # calculating variables
            if (mouse.getPressed()[0]==1) and (kfix == True):
            # if sst: p_port.write(b'31')  # send trigger
                choice_pos=rep_pos
                choice_ang= float("%.2f" % exp.getAngle(rep_pos))
                choice_r=float("%.2f" % (norm(rep_pos)))
                err_ang = np.rad2deg(circdist( np.deg2rad(angle_T),np.deg2rad(choice_ang)))
                err_r=float("%.2f" % (st_prop['radius']-(choice_r)))# -*- coding: utf-8 -*-
            # if sst: p_port.write(b'00')  # send trigger
            else:
                choice_pos=np.array([np.nan,np.nan])
                choice_ang=np.nan  
                choice_r=np.nan
                err_ang=np.nan
                err_r=np.nan

            
            
            if sst: win.callOnFlip(p_port.write, tg_resp.encode()) # send trigger
            #ring.draw()
            routines(False)
            if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # put down pins  
            #ring.draw()
            routines(False)
                            
            
            trial_data = [subj_id, n_trial,thisBlock, t_delay ,kfix, trials[n_trial], last_stim, 
                        deviants, deviant_angles, deviant_pos, deviant_times,  rt, keypressed,
                        r_T,angle_T,resp_mod,
                        choice_pos[0],choice_pos[1],
                        choice_ang,choice_r,
                        resp_time,mtime,
                        m_pos_x,m_pos_y,m_clock]
            if log_on:
                with open(filename, 'a') as csvfile:
                    spamwriter = csv.writer(csvfile, delimiter=';',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    
                    spamwriter.writerow(trial_data)

            # Feedback color scale

            abs_err_angle = np.abs(np.round(err_ang)) 
            if abs_err_angle < 5:
                col_fb = 'green'
            if abs_err_angle >= 5 and abs_err_angle < 25:
                col_fb = 'orange'
            if abs_err_angle >= 25 and abs_err_angle < 60:
                col_fb = 'red'
            if abs_err_angle > 60:
                col_fb = 'black'
            
            if tobii:
                eyes_stamps.append(tr.get_system_time_stamp()) # feedback init time stamp

            for frameN in range(st_prop['FEEDBACK_frames']):
                if kfix == True:
                  #  ring.draw()
                    resp.draw()
                    if ct == True:
                        if keypressed == 'space':
                            msgTex=visual.TextStim(win=win, ori=0,
                            text= 'Detected',
                            pos=[0,-2], height=1,
                            color='green',units="cm")
                            msgTex.draw()
                        else:
                            msgText=visual.TextStim(win=win, ori=0,
                            text= 'Not detected',
                            pos=[0,-2], height=1,
                            color='red',units="cm").draw()
                    if ct == False:
                        if keypressed == 'space':
                            msgText=visual.TextStim(win=win, ori=0,
                            text= 'Detection error',
                            pos=[0,-2], height=1,
                            color='green',units="cm").draw()
                        else:
                            msgText=visual.TextStim(win=win, ori=0,
                            text= abs_err_angle,
                            pos=[0,-2], height=1.5,
                            color=col_fb,units="cm").draw()
                
                else:
                    msgText=visual.TextStim(win=win, ori=0,
                    text= 'KEEP YOUR EYES ON THE FIXATION POINT',
                    pos=[0,-2], height=1,
                    color='red',units="cm").draw()

                routines(False)
            if kfix == True:
                trep = False # end trial and begin the next one

        #exp.say_msg(abs_err_angle,0.3, win)

        #elist[n_trial] = trial_data # append new row data to list
        thr_trials_var.append(trial_data)
        if tobii:
            trial_eYe_data['all_gaze_data'] = alls_gaze_data 
            trial_eYe_data['time_stamps_events'] = eyes_stamps
            eyedata[n_trial] = trial_eYe_data # storing eyedata for this trial

        sleep(np.random.uniform(0.25,0.75)) # add some jitter time between trials

        #resp.setPos(pos)
        #resp.setFillColor('white')
        #resp.draw()
    
    block['data'] = thr_trials_var
    main_exp['Exp_blocks'][thisBlock] =  block
    

    expInfo['main_exp'] =  main_exp

    toFile(filepy, expInfo) #saving file to disk

    if tobii:
        main_eye['Exp_blocks'][thisBlock]  =  eyedata
        toFile(filepeye, expInfo) #saving file to disk
    
            #trialClocktimes.append(trial_times)
        # Get datafiles in pandas format and attack to main Exp.variabl
 
   # if sst: win.callOnFlip(p_port.write, tg_zero.encode()) # send trigger
    routines(False)

if tobii:
    my_eyetracker.unsubscribe_from(tr.EYETRACKER_GAZE_DATA, gaze_data_callback)

main_exp['monitor_features'] = monitor_features
main_exp['stim_props'] = st_prop
#main_exp['lotery'] = corr_lotery

expInfo['main_exp'] =  main_exp

print('Overall, %i frames were dropped.' % win.nDroppedFrames)

toFile(filepy, expInfo) #saving variables to disk    

instr.end_experiment(win)   
#clean up
event.clearEvents()
win.close()                                                              # eyetracker
core.quit
#exp.exit_task(win)



#record onset time for first stim in the sequence
#Update instructions
#write characteristics of the task
#adapt stimuli size and times
#ensure that data are recorded
