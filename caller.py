import PySimpleGUI as sg
from tree_builder import TreeMaker as T
from childnodes import Childnodes as C
'''
layout2 = [[sg.Text('Type YES to produce a phylogenetic tree of entities.')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window2 = sg.Window('Data Entry', layout2)
event, values2 = window2.read()
window2.close()
layout3 = [[sg.Text('Type YES to completely update the existing node descendent file - childnodes.txt (if answered YES previously, else click OK)')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window3 = sg.Window('Data Entry', layout3)
event, values3 = window3.read()
window3.close()
cnodesdir = sg.popup_get_file("Location of existing childnodes.txt file")
sg.popup('You entered', cnodesdir)
layout4 = [[sg.Text('Enter the NCBI REST API key associated with your account')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window4 = sg.Window('Data Entry', layout4)
event, values4 = window4.read()
window4.close()

if type(values2) == dict:
                 
                values2 = str(values2[0])
if type(values3) == dict:
                 
                values3 =  str(values3[0])

if type(values4) == dict:
                 
                values4 = str(values4[0])

if values2 == 'YES':
        ids  = sg.popup_get_file("Location of ids.txt file")
        sg.popup('You entered', ids)
        if values3 == 'YES':
                print("updating the existing childnodes.txt file")
                result2 = C(cnodesdir, ids, values4)
                result2.filerefresh()
        print('Checking if any new species need to be added')
        result2 = C(cnodesdir, ids, values4)
        result2.updatenewspec()
        
        layout5 = [[sg.Text('Enter the tree title')],
        [sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
        window5 = sg.Window('Data Entry', layout5)
        event, values5 = window5.read()
        window5.close()
        outputdir = sg.popup_get_file("Location of existing childnodes.txt file")
        sg.popup('You entered', cnodesdir)
        if type(values5) == dict:
                 
                values5 = str(values5[0])
        result3 = T(ids, outputdir, cnodesdir, values5)
        print('Check complete. Constructing trees.')
        result3.Maker()
'''
#C(childnodesdir, ids found text file, ncbiapikey)

#T(ids found text file, output directory, childnodesdir, title)
result3 = T('Testset4/ids2.txt', 'results', 'childnodes.txt', 'Testset4')
result3.Maker()
