from ete3 import TreeStyle, RectFace, faces, Tree, NodeStyle, TextFace
from ete3 import NCBITaxa
import ast
from Bio import Phylo
import os
import math
# Make sure to install svglib and rlPyCairo - rlPyCairo is not in the imports but is necessary to make svglib and reportlab work
from PIL import Image, ImageFont, ImageDraw
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM
import time
class TreeMaker:
        def __init__(self, items_to_find, directorypath, cnodesdir, treetitle):
            #Initialise inputs
            with open(items_to_find) as sf:
                    text = sf.readlines()
                    newls = []
                    for i in text:
                            newls.append(i.rstrip("\n"))
            self.items_to_find = newls
            self.directorypath = directorypath
            self.cnodesdir = cnodesdir
            self.treetitle = treetitle
        def listmaker(self, listtobeprocessed, allitems):
                #Lists all the descendants of items to be found (the txt file which is the user input) and adds them to allitems
                ncbi = NCBITaxa()
                for itf in listtobeprocessed:
                        for i in ncbi.get_lineage(itf):
                                allitems.append(str(i))
                allitems = list(dict.fromkeys(allitems))
                return allitems
        def colourselecter(self, colourdict):
                #First find all the descendants and taxa of relevance (see listmaker above)
                allitems = self.listmaker(self.items_to_find, [])
                weighteddict=dict.fromkeys(allitems,0)
                ncbi = NCBITaxa()
                #cnodesdir is a txt files of all nodes and how many direct descendants they have. This turns the txt file into a dictionary. 
                with open (self.cnodesdir ,encoding = 'utf-8') as cn:
                        text = cn.readlines()
                        childnodedict = dict()
                        for i in text:
                                item = i.split(" ")
                                if item[1] == "\n":
                                        item[1] = 1
                                else:
                                        item[1] = item[1].strip("\n")
                                key = item[0]
                                value = item[1]
                                childnodedict[key] = value
                 # Iterate through items to find and count.
                itfset = set(self.items_to_find)
                for itf in itfset:
                    weighteddict[itf] = 0
                for itf in self.items_to_find:
                    weighteddict[itf] += 1    
                for key, value in weighteddict.items():
                        if value > 0:
                            weighteddict[key] = math.log10(value)  
                        else:
                            weighteddict[key] = value       
                total = max(weighteddict.values())
                
                viridis = ['#fde725',
'#f8e621',
'#f1e51d',
'#ece51b',
'#e5e419',
'#dfe318',
'#d8e219',
'#d0e11c',
'#cae11f',
'#c2df23',
'#bddf26',
'#b5de2b',
'#addc30',
'#a8db34',
'#a0da39',
'#9bd93c',
'#93d741',
'#8ed645',
'#86d549',
'#7fd34e',
'#7ad151',
'#73d056',
'#6ece58',
'#67cc5c',
'#60ca60',
'#5cc863',
'#56c667',
'#52c569',
'#4cc26c',
'#48c16e',
'#42be71',
'#3dbc74',
'#3aba76',
'#35b779',
'#32b67a',
'#2eb37c',
'#2ab07f',
'#28ae80',
'#25ac82',
'#24aa83',
'#22a785',
'#20a486',
'#1fa287',
'#1fa088',
'#1f9e89',
'#1e9b8a',
'#1f998a',
'#1f968b',
'#20938c',
'#20928c',
'#218f8d',
'#228d8d',
'#238a8d',
'#24878e',
'#25858e',
'#26828e',
'#26818e',
'#277e8e',
'#287c8e',
'#29798e',
'#2a768e',
'#2b748e',
'#2c718e',
'#2d708e',
'#2e6d8e',
'#306a8e',
'#31688e',
'#32658e',
'#33638d',
'#34608d',
'#365d8d',
'#375b8d',
'#38588c',
'#39558c',
'#3b528b',
'#3c508b',
'#3d4d8a',
'#3e4989',
'#3f4788',
'#414487',
'#424186',
'#433e85',
'#443a83',
'#453882',
'#463480',
'#46327e',
'#472e7c',
'#472c7a',
'#482878',
'#482475',
'#482173',
'#481d6f',
'#481b6d',
'#481769',
'#471365',
'#471063',
'#460b5e',
'#46085c',
'#450457',
'#440154']
                reverseviridis = viridis[::-1]
                #Attributes a number to each colour in the viridis scale for accessing later. 
                colourscaledict = {}
                for i in range(0, len(reverseviridis)):
                        colourscaledict[str(reverseviridis[i])]= i
                tblabelled = []
                # Deciding which taxa should be annotated with names (and added to the tblabelled list)
                for key in weighteddict.keys():
                        rankdict = ncbi.get_rank([key]) 
                        nottblabelled = False
                        if 'subspecies' in rankdict.values() and weighteddict[key]*100/total>10:
                                lineage = ncbi.get_lineage(key)
                                LEN = len(lineage)
                                parent = lineage[LEN-1]
                                parentweight = weighteddict[str(parent)]
                                descendants = ncbi.get_descendant_taxa(key)
                                descendantweight = 0
                                for i in descendants:
                                        if i in weighteddict.keys():
                                                descendantweight += weighteddict[str(i)]
                                        else:
                                                continue
                                if weighteddict[key]>parentweight and weighteddict[key]>descendantweight:
                                        tblabelled.append(key)
                        elif 'species' in rankdict.values() and weighteddict[key]*100/total >40 :
                                lineage = ncbi.get_lineage(key)
                                LEN = len(lineage)
                                parent = lineage[LEN-1]
                                parentweight = weighteddict[str(parent)]
                                descendants = ncbi.get_descendant_taxa(key)
                                descendantweight = 0
                                for i in descendants:
                                        if i in weighteddict.keys():
                                                descendantweight += weighteddict[str(i)]
                                        else:
                                                continue
                                if weighteddict[key]>parentweight and weighteddict[key]>descendantweight:
                                        tblabelled.append(key)
                                
                        elif 'genus' in rankdict.values() and weighteddict[key]*100/total>50:
                                lineage = ncbi.get_lineage(key)
                                LEN = len(lineage)
                                parent = lineage[LEN-1]
                                parentweight = weighteddict[str(parent)]
                                descendants = ncbi.get_descendant_taxa(key)
                                descendantweight = 0
                                for i in descendants:
                                        if i in weighteddict.keys():
                                                descendantweight += weighteddict[str(i)]
                                        else:
                                                continue
                                if weighteddict[key]>parentweight or weighteddict[key]>descendantweight:
                                        tblabelled.append(key)
                        elif 'family' in rankdict.values() and weighteddict[key]*100/total>60:
                                lineage = ncbi.get_lineage(key)
                                LEN = len(lineage)
                                parent = lineage[LEN-1]
                                parentweight = weighteddict[str(parent)]
                                descendants = ncbi.get_descendant_taxa(key)
                                descendantweight = 0
                                for i in descendants:
                                        if i in weighteddict.keys():
                                                descendantweight += weighteddict[str(i)]
                                        else:
                                                continue
                                if weighteddict[key]>parentweight and weighteddict[key]>descendantweight:
                                        tblabelled.append(key)
                        
                          
                        
                for key,value in weighteddict.items():
                        if value == 0:
                                colourdict[key] = '#440154'
                        else:       
                                placeindex = (value / total) * 100
                                placeindex = math.ceil(placeindex)
                                colourdict[key] = reverseviridis[placeindex - 1]
                penultimate = []
                for tbl in tblabelled:
                        lineagetbl = ncbi.get_lineage(tbl)
                        penultimate.append(lineagetbl[len(lineagetbl)-2])
                parents = set(penultimate)
                parents = list(parents)
                siblings = {}
                for p in parents:
                        sib = []
                        for tbl in tblabelled:
                                lineagetbl = ncbi.get_lineage(tbl)
                                parent = (lineagetbl[len(lineagetbl)-2])
                                if parent == p:
                                        sib.append(tbl)
                        siblings[p] = sib
                for tbl in tblabelled:
                        lineagetbl = ncbi.get_lineage(tbl)
                        parent = (lineagetbl[len(lineagetbl)-2])
                        sibs = siblings[parent]
                        heaviest = 0
                        if len(sibs) >1:
                                sibweights = []
                                for s in sibs:
                                        sibweights.append(weighteddict[s])
                                heaviest = max(sibweights)
                        if weighteddict[tbl] < heaviest:
                                tblabelled.remove(tbl)



                #Save important lists and dictionaries to txtfiles folder
                os.mkdir(self.directorypath + "/txtfiles")
                time.sleep(2)
                with open(self.directorypath + '/txtfiles/colourdict.txt', 'w', encoding = 'utf-8') as f:
                        f.write(str(colourdict))
                with open(self.directorypath + '/txtfiles/total.txt', 'w', encoding = 'utf-8') as g:
                        g.write(str(total))
                with open(self.directorypath + '/txtfiles/tblabelled.txt', 'w', encoding = 'utf-8') as h:
                        h.write(str(tblabelled))
                with open(self.directorypath + '/txtfiles/weighteddict.txt', 'w', encoding = 'utf-8') as i:
                        i.write(str(weighteddict))
                #Pause function to prevent the errno 13 error
                time.sleep(2)
                with open(self.directorypath + '/txtfiles/colourscaledict.txt', 'w', encoding = 'utf-8') as j:
                        j.write(str(colourscaledict))
                time.sleep(2)
                
                                        

        def layoutfunc(self, node):
                  ncbi = NCBITaxa()
                  rankdict = ncbi.get_rank([node.name])

                  node.complete_branch_lines_when_necessary = False
                  node.optimal_scale_level = "full"
                  node.guiding_lines_type = 0
                  node.extra_branch_line_type = 0
                  with open(self.directorypath + '/txtfiles/colourdict.txt') as f:
                                colourdict = f.read()
                  with open(self.directorypath + '/txtfiles/tblabelled.txt') as g:
                                tblabelled = g.read()
                  with open(self.directorypath + '/txtfiles/colourscaledict.txt', 'r') as j:
                                indicator = j.read()             
                  colourdict = ast.literal_eval(colourdict)
                  tblabelled = ast.literal_eval(tblabelled)
                  indicator = ast.literal_eval(indicator)
                  tblabellednames = ncbi.get_taxid_translator(tblabelled)
                  
                  
                  node.img_style["hz_line_type"] = 0
                  
                  if node.get_children() == [] or node.name not in colourdict.keys():
                          node.img_style["hz_line_color"] = "#ffffff"
                  nohorline = False
                  

                  if node.name in colourdict.keys():
                        
                        if indicator[colourdict[str(node.name)]] >= 60:
                                amplifier = 15
                        elif indicator[colourdict[str(node.name)]] <=59 and indicator[colourdict[str(node.name)]] >=20:
                                amplifier = 10
                        elif indicator[colourdict[str(node.name)]] <=19 and indicator[colourdict[str(node.name)]] >=3:
                                amplifier = 3
                        else:
                                amplifier = indicator[colourdict[str(node.name)]]
                        node.img_style["fgcolor"] = colourdict[str(node.name)]
                        node.img_style["size"] =0
                        node.img_style["vt_line_color"] = colourdict[str(node.name)]
                        node.img_style["hz_line_color"] = colourdict[str(node.name)]
                        node.img_style["vt_line_width"] = 0
                        node.img_style["hz_line_width"] = 0

                        if node.get_children == []:
                                for i in node.get_children():
                                        if i in colourdict.keys():
                                                continue
                                        else:
                                                nohorline == True
                        
                        if node.name in tblabelled:
                                        if 'species' or 'subspecies' or 'genus' in rankdict.values():
                                                        node.img_style["hz_line_type"] = 0
                                        if amplifier == 15:
                                                textsize = 10
                                                node.img_style["hz_line_type"] = 1
                                                faces.add_face_to_node(TextFace(tblabellednames[int(node.name)], ftype = 'arial', fsize = textsize, fgcolor = "000000", penwidth=0, fstyle= 'normal', tight_text = False, bold = False), node, column = 1, position = "float")
                                        elif amplifier == 10:
                                                textsize = 9
                                                node.img_style["hz_line_type"] = 1
                                                faces.add_face_to_node(TextFace(tblabellednames[int(node.name)], ftype = 'arial', fsize = textsize, fgcolor = "000000", penwidth=0, fstyle= 'normal', tight_text = False, bold = False), node, column = 3, position = "float")
                                        elif amplifier == 3:
                                                textsize = 8
                                                node.img_style["hz_line_type"] = 1
                                                faces.add_face_to_node(TextFace(tblabellednames[int(node.name)], ftype = 'arial', fsize = textsize, fgcolor = "000000", penwidth=0, fstyle= 'italic', tight_text = False, bold = False), node, column = 3, position = "float")
                                        else:
                                                textsize = 0 
                  else: 
                        
                        node.img_style["size"] = 0
                        node.img_style["fgcolor"] = "#ffffff"
                        node.img_style["bgcolor"] = "#ffffff"
                        node.support = 0
                        node.distance = 0
                        node.img_style["vt_line_color"] = "#ffffff"
                        node.img_style["hz_line_color"] = "#ffffff"
                        
                        if node.get_children() != []:
                                nohorline == True
                                        
                  if nohorline == True:
                              node.img_style["hz_line_color"] = "#ffffff"   
                            
        def Maker(self):     
                                
                                     with open(self.cnodesdir,encoding = 'utf-8') as fp:
                                                      
                                                      text = fp.readlines()
                                                      topologyfeeder = []
                                                      for i in text:
                                                            if " " in i:
                                                                    
                                                                id = i.split(" ")[0]
                                                            else:
                                                                id = i
                                                            topologyfeeder.append(str(id))
                                                
                                                      ncbi = NCBITaxa()
                                                      
                                                      tree = ncbi.get_topology(topologyfeeder, intermediate_nodes=True)
                                                      tree.annotate_ncbi_taxa()
                                                      ts = TreeStyle()
                                                      ts.layout_fn = self.layoutfunc
                                                      ts.show_leaf_name = False
                                                      ts.mode = "c"
                                                      ts.root_opening_factor = 0
                                                      ts.arc_start = 0 # 0 degrees = 3 o'clock
                                                      ts.arc_span = 360
                                                      self.colourselecter({})
                                                      tree.show(tree_style=ts)
                                                      os.mkdir(self.directorypath + "/trees")
                                                      tree.write(format = 0, outfile = self.directorypath + "/trees/new_tree.nwk")
                                                      tree.render(self.directorypath + "/trees/img_tree.svg", w= 3600, units = 'px', tree_style = ts)
                                                      Phylo.convert(self.directorypath + "/trees/new_tree.nwk", "newick", self.directorypath + "/trees/new_tree.xml", "nexml")
                                                      img = svg2rlg(self.directorypath + "/trees/img_tree.svg")
                                                      renderPM.drawToFile(img, self.directorypath + "/trees/tree.png", fmt = "PNG")

                                                      filename = self.directorypath + "/trees/tree.png"
                                                      with Image.open(filename) as img: 
                                                                width, height = img.size
                                                                img = img.resize((width * 2, height * 2 ))
                                                                img.save(self.directorypath + "/trees/tree.png")
                                                                img2 = Image.open('colourbar.png')
                                                                img.paste(img2, (10, 10))
                                                                img.save(self.directorypath + "/trees/tree+bar.png")
                                                      with open(self.directorypath + '/txtfiles/total.txt') as j:
                                                                        total = j.read()
                                                      img = Image.open(self.directorypath + "/trees/tree+bar.png")
                                                      draw = ImageDraw.Draw(img)
                                                      font = ImageFont.truetype("arial", 50)
                                                      font2 = ImageFont.truetype("arial", 70)
                                                      draw.text((600, 150), "Weighted total=  " + str(total), (0, 0, 0), font = font)
                                                      
                                                      img.save(self.directorypath + "/trees/tree+bar.png")
                                                      draw.text((2500, 150), "Phylogenetic Tree based on the named entities in the folder: "  + self.treetitle, (0, 0, 0), font = font2)
                                                      img.save(self.directorypath + "/trees/tree+bar.png")