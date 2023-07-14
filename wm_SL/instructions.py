

from psychopy import visual, event, core
import random
import numpy as np
# Experiment instructions




def main_instructions(win):
       
    inst = visual.TextStim(win, pos = [0,2])
   
    inst.text = "Practica:\
    En este experimento, mientras mantienes tus ojos en el punto de fijación central, tendrás que prestar atención a una serie de circulos que apareceran rápidamente alrededor.\
    Tu tarear consistirá en:\
    1) pulsar la barra espaciadora con la mano izquierda tan pronto como detectes que alguno de los estímulos es diferente al resto.\
    2) al final de la secuencia, arrastrar el ratón hasta la posición en la que viste el último estimulo y hacer click."
    inst.height = 0.7

    nextt = visual.TextStim(win, pos = [0,-8])
    nextt.wrapWidth = 20
    nextt.height = 0.7
    nextt.color = 'black'
    nextt.text = "Press space to continue "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])

    
    inst.text = "Practica:\
    Trás señalar la posición del último estímulo, aparecerá un número indicando cuan lejos estaba tu respuesta de la posición real del estímulo. \
    Cuanto menor sea esta diferencia mejor.\
    Recuerda que durante cada trial, deberás mantener la mirada en el punto de fijación (y puedes descansar al final de cada trial o bloque) \
    El experimento tiene una duración de 8 bloques. En total son 60 minutos. Muchas gracias por participar!!"
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    
    return

def block_start(win, nblock, blocklen, cond_resps):
    inst = visual.TextStim(win, pos = [0,-4])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'black'
    nextt = visual.TextStim(win, pos = [0,-8])
    nextt.height = 0.7
    nextt.color = "black"
    
    mod = visual.TextStim(win, pos = [0, 5])
    mod.height = 1.2
    mod.color = "white"
    
    inst.text = 'Comenzando bloque ' + str(nblock+1) + ' / ' + str(blocklen)
    nextt.text = "Cuando estes list@ pulsa la barra espaciadora "
    inst.draw()
    nextt.draw()
    win.flip()
    event.waitKeys( keyList=['space'])
    core.wait(1)
    return

    
def new_trial(win):
    inst = visual.TextStim(win, pos = [0,0])
    inst.wrapWidth = 20
    inst.height = 0.7
    inst.color = 'white'
    inst.text = "New trial"
    inst.draw()
    win.flip()
    core.wait(0.75)
    return

    


def end_experiment(win):
    inst = visual.TextStim(win, pos=[0,0], height = 1.2)
    inst.text = "Final del experimento! Avisa al investigador"     
    inst.draw()
    win.flip()
    event.waitKeys(keyList = ["space"])
    return

        