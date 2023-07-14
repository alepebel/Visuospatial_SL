from psychopy import visual
import numpy as np
#import orient_decision.win
 
def stim_config(ifi):
    # Experiment timings design
    
    stim = {}
    stim['CUE_time'] = 150 # in miliseconds
    stim['CUE_frames'] = round(stim['CUE_time']/ifi)
    stim['FP'] = round(750/ifi) # fixation  + circle + resp_options
    stim['FEEDBACK_frames'] = round(500/ifi)
    #stim['ENDTRIAL_frames'] = round(500/ifi) # fixation  + circle + resp_options

    # Create all the stimuli
    # Stimuli parameters
    stim['size_CUE'] = 2 # degs # este es el lado del cuadrado del grating
    stim['radius'] = 10

    return stim

def draw_basic(win, stim): # creates a dictionary with basic stimulus features
    
    # Fication point
    basic_stim = {}
    basic_stim['fixation_point'] =  visual.PatchStim(win, color= [1, 1, 1], tex=None,mask='circle', size=0.3, units = "deg")
     
    # Grating stim
    basic_stim['grating'] = visual.GratingStim(win=win, mask="raisedCos" , size=stim['size_stim'], pos=[0,0], sf=stim['SF'],
                                 units = "deg", contrast = stim['grating_contrast'], maskParams = {'sd': 3} ) # , color = [1,0,1]
    
    # Contour boundary
    basic_stim['stim_contour'] = visual.Circle(win=win,lineWidth = 10, units="deg", radius=stim['size_stim']/2, lineColor=[-1, -1, -1],edges=128)
    basic_stim['stim_contour_in'] = visual.Circle(win=win,lineWidth = 4, units="deg", radius=stim['size_stim']/2, lineColor=[1, 1, 1],edges=128)

    return basic_stim


def fixation(win, basic_stim):
    basic_stim['fixation_circle'].draw()
    basic_stim['fixation_point_c'].draw()
    basic_stim['fixation_point'].draw()
    return

def draw_contour(win, basic_stim):
    basic_stim['stim_contour'].draw()
    basic_stim['stim_contour_in'].draw()
    return



def monitor_def():
     # lets define the monitor if necessary in a dictionary (you can define other monitors def2, def3)
    monitor_def1 = {}
    monitor_def1['monitor_name'] = 'mundet_screen' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def1['monitor_pixels']  = (1920, 1080)
    monitor_def1['monitor_width'] = 53
    monitor_def1['distance2monitor'] = 50
    
    monitor_def2 = {}
    monitor_def2['monitor_name'] = 'multiple1' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def2['monitor_pixels']  = (1920, 1080)
    monitor_def2['monitor_width'] = 53
    monitor_def2['distance2monitor'] = 50
    
    monitor_def3 = {}
    monitor_def3['monitor_name'] = 'asus_home' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def3['monitor_pixels']  = (1920,1080)
    monitor_def3['monitor_width'] = 60
    monitor_def3['distance2monitor'] = 40
    
    monitor_def4 = {}
    monitor_def4['monitor_name'] = 'IDIBAPS' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def4['monitor_pixels']  = (1920,1080)
    monitor_def4['monitor_width'] = 58
    monitor_def4['monitor_height'] = 26.5
    monitor_def4['distance2monitor'] = 60
    
    monitor_def5 = {}
    monitor_def5['monitor_name'] = 'mac_pro13' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def5['monitor_pixels']  = (1200, 800) #(2560,1600)
    monitor_def5['monitor_width'] = 29
    monitor_def5['monitor_height'] = 18
    monitor_def5['distance2monitor'] = 40

    monitor_def6 = {}
    monitor_def6['monitor_name'] = 'LENOVO_IDI' # monitor to use (make sure to define monitors in exp_monitors center.
    monitor_def6['monitor_pixels']  = (1920,1080)
    monitor_def6['monitor_width'] = 38
    monitor_def6['monitor_height'] = 21.15
    monitor_def6['distance2monitor'] = 40
    
    
    monitors = [monitor_def1, monitor_def2, monitor_def3, monitor_def4, monitor_def5, monitor_def6]
    return monitors
        
        