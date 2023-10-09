import PySimpleGUI as sg
import childnodes as C
import Treebuilder as T
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
        print('Check complete. Constructing trees.')
        result3 = T(ids, outputfiles, cnodesdir, inputfiles)
        result3.Maker()
'''
# ann(tax dictionary, input directory, output directory, count, section, casesens)
#result = ann('NCBI_tax_dictionary8.json', 'testset', 'results', 0, 'ALL', 'NO')
#result.initialsteps()
# C(childnodesdir, ids found text file, ncbiapikey)
#result2 = C('childnodes.txt', 'results/ALLAnnotated_output_2023-09-22_17-01-00/ids.txt', 'f55726c2c32772c2b82304814b30148aff07')
#result2.updatenewspec()
# T( ids found text file, output directory, childnodesdir, input directory)
result3 = T('Testset4/ids2.txt', 'results4', 'childnodes.txt', 'Testset4')
result3.Maker()
