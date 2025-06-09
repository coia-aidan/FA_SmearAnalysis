#!/usr/bin/env python3

import PySimpleGUI as sg
import CircScript2 as script


layout = [[sg.Text("Prepare files in ProSize: ",s=(80,1),font='Georgia',text_color='cyan2',background_color='midnight blue')],
	  [sg.Text("Ensure min. peak height is set below peaks (2 intron peaks, main peak, main concatemer peak.",size=(120,1),font='Georgia',text_color='cyan2',background_color='midnight blue')],
	  [sg.Text("Export run data (With Alternate Peak Table selected)",size=(80,2),font='Georgia',text_color='cyan2',background_color='midnight blue')],
	  [sg.HorizontalSeparator(color = 'black',pad = ((10,10),(20,20)))],
	  
	  [sg.FileBrowse("Select Peak Table .csv: ",key='-PEAK-',size = (40, 2),button_color = 'LightSteelBlue4',font = 'Georgia',pad = ((20,20),(1,1)))],
	  [sg.Text(size=(80,2), key='-OUTPUT1-',font = 'Beirut',text_color = 'black', background_color = 'grey92')],
	 
	  [sg.FileBrowse("Select ANAI File: ",key='-ANAI-',size = (40, 2),button_color = 'SkyBlue4',font = 'Georgia',pad = ((20,20),(1,1)))],
	  [sg.Text(size=(80,2), key='-OUTPUT2-',font = 'Beirut',text_color = 'black', background_color = 'grey99')],
          
 	  [sg.HorizontalSeparator(color = 'black',pad = ((10,10),(20,20)))],
	  [sg.Button('Run', button_color=('midnight blue','green3'), border_width=10, s=(20,2),mouseover_colors=('lawn green','lawn green'),font='Beirut'), 
	   sg.Button('Exit', button_color=('firebrick4','tomato2'), border_width=10, s=(20,2),mouseover_colors=('firebrick1','firebrick1'),font='Beirut')]]


window = sg.Window('CircRNA Smear Analysis Version 2', layout,resizable=True, text_justification = 'center',
			element_justification = "center", titlebar_font = 'Trattatello',
			element_padding = 2, background_color = "DeepSkyBlue3",
			margins = (5, 5), font = "SignPainter") 


while True:
    event, values = window.read()
    if event == sg.WINDOW_CLOSED or event == 'Exit':
        break
    window['-OUTPUT1-'].update(values['-PEAK-'])
    window['-OUTPUT2-'].update(values['-ANAI-'])

    response = script.inputs(values['-PEAK-'], values['-ANAI-']), print(values['-PEAK-'], values['-ANAI-'])
    sg.PopupAutoClose("Success. ANAI Updated", no_titlebar=True, background_color = None, 
			text_color = 'green', auto_close = True, auto_close_duration = 30)

window.close()








# Define the window's contents
# Create the window
# [sg.Text("Peak Table .csv: ")],
  #        [sg.Input(key='-INPUT1-')],
   #       [sg.Text("ANAI file .anai: ")],
    #      [sg.Input(key='-INPUT2-')],
# Display and interact with the Window using an Event Loop
  # See if user wants to quit or window was closed
 

# run FAscript.py with inputs - values['-INPUT1-'], values['-INPUT2-']

# if script returns True -> popup success message
    
    # Output a message to the window

#IF successful smear analysis:

   
# Finish up by removing from the screen




