from ete3 import TreeStyle, RectFace, faces, Tree, NodeStyle
from ete3 import NCBITaxa

from Bio import Phylo
import os
import math
# Make sure to install svglib and rlPyCairo - rlPyCairo is not in the imports but is necessary to make svglib and reportlab work
from PIL import Image, ImageFont, ImageDraw
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM

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
                ncbi = NCBITaxa()
                for itf in listtobeprocessed:
                        for i in ncbi.get_lineage(itf):
                                allitems.append(str(i))
                allitems = list(dict.fromkeys(allitems))
                return allitems
        def colourselecter(self, colourdict):
                allitems = self.listmaker(self.items_to_find, [])
                iddict=dict.fromkeys(allitems,0)
                ncbi = NCBITaxa()
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
                for itf in self.items_to_find:
                        reversedls = ncbi.get_lineage(itf)[::-1]
                        factor = 1
                        for index, i in enumerate(reversedls):
                                        if index == 0:
                                                iddict[str(i)] += 1 
                                        elif str(i) in childnodedict.keys():

                                                        newfactor = int(childnodedict[str(i)])
                                                        factor = factor / newfactor
                                                        iddict[str(i)] += factor
                                        else:
                                                        iddict[str(i)] += factor 
                total = max(iddict.values())
                
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
                heavy = []
                for key in iddict.keys():
                        rankdict = ncbi.get_rank([key]) 
                        if 'species' in rankdict.values() and iddict[key] > 2.5 :
                                length = len(ncbi.get_lineage(key))
                                parent = ncbi.get_lineage(key)[length - 2]
                                if parent in heavy:
                                        heavy.remove(parent)
                                        heavy.append(key)
                        elif 'subspecies' in rankdict.values() and iddict[key]>2:
                                heavy.append(key)
                        elif 'genus' in rankdict.values() and iddict[key]>4:
                                heavy.append(key)
                        elif 'family' in rankdict.values() and iddict[key]>5:
                                heavy.append(key)
                        elif 'kingdom' in rankdict.values():
                                heavy.append(key)
                        elif 'domain' in rankdict.values():
                                heavy.append(key)
                        elif iddict[key]> 0.7*total:
                                heavy.append(key)       
                        
                for key,value in iddict.items():
                        placeindex = (value / total) * 100
                        placeindex = math.ceil(placeindex)
                        colourdict[key] = reverseviridis[placeindex - 1]

                delivery = [colourdict, total, heavy, iddict]
                return(delivery) 
                                        

        def layoutfunc(self, node):
                  ncbi = NCBITaxa()
                  node.complete_branch_lines_when_necessary = False
                  node.optimal_scale_level = "full"
                  node.guiding_lines_type = 0
                  node.extra_branch_line_type = 0
                  delivery = self.colourselecter({})
                  colourdict = delivery[0]
                  tblabelled = delivery[2]
                  tblabellednames = ncbi.get_taxid_translator(tblabelled)
                  node.img_style["hz_line_type"] = 0
                  if node.get_children() == [] or node.name not in colourdict.keys():
                          node.img_style["hz_line_color"] = "#ffffff"
                  nohorline = False
                  
                  
                  if node.name in colourdict.keys():
                        node.img_style["fgcolor"] = colourdict[str(node.name)]
                        node.img_style["vt_line_color"] = colourdict[str(node.name)]
                        node.img_style["hz_line_color"] = colourdict[str(node.name)]
                        node.img_style["vt_line_width"] = 5
                        node.img_style["hz_line_width"] = 5
                        if node.get_children == []:
                                for i in node.get_children():
                                        if i in colourdict.keys():
                                                continue
                                        else:
                                                nohorline == True
                        if node.name in tblabelled:
                                #print(tblabellednames[int(node.name)])
                                faces.add_face_to_node(RectFace(0.1, 0.1, fgcolor = "000000", bgcolor= "000000", label = {'text': tblabellednames[int(node.name)] , 'font': 'arial', 'fontsize' : 25}), node, column = 1, position = "float")      
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
                                                            id = i.split(" ")[0]
                                                            
                                                            topologyfeeder.append(str(id))
                                                      #print(topologyfeeder)
                                                
                                                      ncbi = NCBITaxa()
                                                      
                                                      tree = ncbi.get_topology(topologyfeeder, intermediate_nodes=True)
                                                      tree.annotate_ncbi_taxa()
                                                      #print(tree.get_ascii(attributes=["sci_name", "rank"]))
                                                      ts = TreeStyle()
                                                     
                                                      ts.show_leaf_name = False
                                                      ts.mode = "c"
                                                      ts.root_opening_factor = 0
                                                      ts.arc_start = -180 # 0 degrees = 3 o'clock
                                                      ts.arc_span = 360
                                                      ts.layout_fn = self.layoutfunc
                                                      alltaxaintext = self.colourselecter({})[3]
                                                      table = ""
                                                      for key, value in alltaxaintext.items():
                                                              val = str(round(value, 2))
                                                              table += "Taxa " + key + " " + "Weighting " + val + "\n"
                                                      print(table)
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
                                                      total = self.colourselecter({})[1]
                                                      img = Image.open(self.directorypath + "/trees/tree+bar.png")
                                                      draw = ImageDraw.Draw(img)
                                                      font = ImageFont.truetype("arial", 50)
                                                      font2 = ImageFont.truetype("arial", 70)
                                                      draw.text((600, 150), "Weighted total=  " + str(total), (0, 0, 0), font = font)
                                                      
                                                      img.save(self.directorypath + "/trees/tree+bar.png")
                                                      draw.text((2500, 150), "Phylogenetic Tree based on the named entities in the folder: "  + self.treetitle, (0, 0, 0), font = font2)
                                                      img.save(self.directorypath + "/trees/tree+bar.png")