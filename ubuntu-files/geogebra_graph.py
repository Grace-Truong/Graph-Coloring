# Sage code for converting a GeoGebra figure into a SageMath graph, for plotting the resulting SageMath graph, and for exporting a SageMath graph into a GeoGebra figure.

import zipfile
import xml.etree.ElementTree as ET

import sage.all


def geogebra_to_graph(ggb_filename):
    '''
    Converts a GeoGebra file into a SageMath graph.
    The following conventions are followed:
    The GeoGebra figure should be a 2-dimensional geometry construction.
    Points become vertices.
    Line segments become edges [be careful that *lines* are not interpreted as edges, and line segments that are drawn over a point does not connect to that point.]
    
    Various attributes are also stored:
    For points/vertices: label, color, size, style.
    For line segments/edges: label, color, thickness.
    
    The position of the points is also stored in the graph's pos dictionary for later plotting.
    '''
    
    with zipfile.ZipFile(ggb_filename,'r') as f:
        root=ET.fromstring(f.read("geogebra.xml"))
        # the root tag is geogebra
    
    c=root.find('construction')

    vertex_labels=[]  # vertices, as labels
    vertex_information=[]  # each vertex will have a dict of information (such as label, color, etc) about that vertex
    v=0  # current vertex number
    pos=dict()

    # interpret points as vertices
    for elt in c.findall('element'):
        if elt.attrib['type']=='point':
            vertex_labels.append(elt.attrib['label'])
            coords=elt.find('coords')
            pos[v]=(float(coords.attrib['x']),float(coords.attrib['y']))
            
            # save vertex information, such as label, color, etc, in a dict
            info=dict()
            info['label']=vertex_labels[-1]
            color=elt.find('objColor')
            info['color']=(int(color.attrib['r']),
                           int(color.attrib['g']),
                           int(color.attrib['b']),
                           float(color.attrib['alpha']))
            info['size' ]=int(elt.find('pointSize').attrib['val'])
            info['style']=int(elt.find('pointStyle').attrib['val'])
            vertex_information.append(info)
            
            v+=1

    G=sage.graphs.graph.Graph()
    G.add_vertices(range(v))
    G.set_pos(pos)
    
    # store information for each vertex in the Sage graph
    for i in range(v):
        G.set_vertex(i,vertex_information[i])

    # interpret line segments as edges
    for elt in c.findall('command'):
        if elt.attrib['name']=='Segment':
            e=[]
            endpts=elt.find('input')
            e.append(vertex_labels.index(endpts.attrib['a0']))  # first endpoint of this edge
            e.append(vertex_labels.index(endpts.attrib['a1']))  # second endpoint of the edge
            label=elt.find('output').attrib['a0']
            G.add_edge(e)
            G.set_edge_label(e[0],e[1],label)
    
    # we need to extract the edge properties from a different object
    for e in G.edges(sort=False):
        label=G.edge_label(e[0],e[1])
        
        for elt in c.findall('element'):
            if elt.attrib['type']=='segment' and elt.attrib['label']==label:
                # save vertex information, such as label, color, etc, in a dict
                info=dict()
                info['label']=label
                color=elt.find('objColor')
                info['color']=(int(color.attrib['r']),
                               int(color.attrib['g']),
                               int(color.attrib['b']),
                               float(color.attrib['alpha']))
                info['thickness']=int(elt.find('lineStyle').attrib['thickness'])
                G.set_edge_label(e[0],e[1],info)
    
    return G


def geogebra_graph_plot(G):
    '''Returns the plot of graph G, using the information retrieved from the GeoGebra file and stored in G.'''
    
    colors=dict()
    # plotted colors are named colors in matplotlib
    # https://matplotlib.org/3.1.0/gallery/color/named_colors.html
    colors[(0,255,0)]='springgreen'  # green in geogebra; vertices to extend coloring to
    colors[(255,0,0)]='salmon'  # red in geogebra; vertex and edge reducers
    colors[(0,255,255)]=other_v='skyblue'  # the rest of the vertices
    colors[(0,0,0)]=other_e='black'  # the rest of the edges
    colors[(255,0,255)]='magenta'  # magenta in geogebra; vertices to recolor
    
    vertex_colors=dict()
    for color in colors.values():
        vertex_colors[color]=[]
    for v in G.vertices(sort=False):
        try:
            color_of_v=G.get_vertex(v)['color'][:3]
        except TypeError:  # label is not a dictionary, or does not have 'color'
            color_of_v=other_v
        if color_of_v in colors:
            vertex_colors[colors[color_of_v]].append(v)
        else:
            vertex_colors[other_v].append(v)
    
    edge_colors=dict()
    for color in colors.values():
        edge_colors[color]=[]
    for e in G.edges(sort=False):
        try:
            color_of_e=G.edge_label(e[0],e[1])['color'][:3]
        except TypeError:  # label is not a dictionary, or does not have 'color'
            color_of_e=other_e
        if color_of_e in colors:
            edge_colors[colors[color_of_e]].append(e)
        else:
            edge_colors[other_e].append(e)
    
    return(G.plot(vertex_colors=vertex_colors,edge_colors=edge_colors))


def graph_to_geogebra(G,ggb_filename):
    '''
    Converts a SageMath graph into a GeoGebra file, where vertices become points and edges become line segments.
    Vertex labels become vertex labels.
    The graph's pos dictionary is used for the position of the points.
    '''
    with open("geogebra.xml","wt") as f:
        f.write(r'<?xml version="1.0" encoding="utf-8"?>' +"\n")
        f.write(r'<geogebra format="5.0">' +"\n")
        
        f.write(r'<euclidianView>' +"\n")
        f.write(r'<viewNumber viewNo="1"/>' +"\n")
        #f.write(r'<size width="2339" height="1196"/>' +"\n")
        #f.write(r'<coordSystem xZero="888.4649186451074" yZero="465.3980334757529" scale="653.3442584680625" yscale="653.3442584680623"/>' +"\n")
        f.write(r'<evSettings axes="false" grid="false" gridIsBold="false" pointCapturing="3" rightAngleStyle="1" checkboxSize="26" gridType="3"/>' +"\n")
        f.write(r'<bgColor r="255" g="255" b="255"/>' +"\n")
        f.write(r'<axesColor r="0" g="0" b="0"/>' +"\n")
        f.write(r'<gridColor r="192" g="192" b="192"/>' +"\n")
        f.write(r'<lineStyle axes="1" grid="0"/>' +"\n")
        f.write(r'<axis id="0" show="true" label="" unitLabel="" tickStyle="1" showNumbers="true"/>' +"\n")
        f.write(r'<axis id="1" show="true" label="" unitLabel="" tickStyle="1" showNumbers="true"/>' +"\n")
        f.write(r'</euclidianView>' +"\n")
        
        f.write(r'<construction>' +"\n")
        
        pos=G.get_pos()
        for v in G.vertices(sort=False):
            # see if there is an object associated with v; if so, interpret it as a dictionary
            info=G.get_vertex(v)
            if info==None:
                info=dict()  # empty dictionary
            
            f.write(f'<element type="point" label="{v}">' +"\n")  # cannot use info['label'], since this needs to match with the Sage vertex name for the edges.
            f.write(r'<show object="true" label="true"/>' +"\n")
            if 'color' in info:
                f.write(f"<objColor r=\"{info['color'][0]}\""
                                +f" g=\"{info['color'][1]}\""
                                +f" b=\"{info['color'][2]}\""
                            +f" alpha=\"{info['color'][3]}\"/>" +"\n")
            else:  # default color
                f.write(r'<objColor r="77" g="77" b="255" alpha="0"/>' +"\n")
            f.write(r'<layer val="0"/>' +"\n")
            f.write(r'<labelMode val="0"/>' +"\n")
            f.write(f'<coords x="{float(pos[v][0])}" y="{float(pos[v][1])}" z="1"/>' +"\n")  # force evaluation of expressions as floats
            if 'size' in info:
                f.write(f"<pointSize val=\"{info['size']}\"/>" +"\n")
            else:
                f.write(r'<pointSize val="5"/>' +"\n")
            if 'style' in info:
                f.write(f"<pointStyle val=\"{info['style']}\"/>" +"\n")
            else:
                f.write(r'<pointStyle val="0"/>' +"\n")
            f.write(r'</element>' +"\n")
        
        for e in G.edges(sort=False):
            label=f"edge-{e[0]}-{e[1]}"
            f.write(r'<command name="Segment">' +"\n")
            f.write(f'<input a0="{e[0]}" a1="{e[1]}"/>' + "\n")
            f.write(f'<output a0="{label}"/>' +"\n")
            f.write(r'</command>' +"\n")
            f.write(f'<element type="segment" label="{label}">' +"\n")
            f.write(r'<show object="true" label="false"/>' +"\n")  # do not show edge label
            f.write(r'<objColor r="0" g="0" b="0" alpha="0"/>' +"\n")
            f.write(r'<layer val="0"/>' +"\n")
            f.write(r'<labelMode val="0"/>' +"\n")
            #f.write(r'<coords x="-0.8400000000000001" y="1.2800000000000002" z="4.5984"/>' +"\n")
            f.write(r'<lineStyle thickness="5" type="0" typeHidden="1" opacity="178"/>' +"\n")
            #f.write(r'<outlyingIntersections val="false"/>' +"\n")
            #f.write(r'<keepTypeOnTransform val="true"/>' +"\n")
            f.write(r'</element>' +"\n")
        
        f.write(r'' +"\n")
        f.write(r'</construction>' +"\n")
        f.write(r'</geogebra>' +"\n")
        # f.write(r'' +"\n")
    
    with zipfile.ZipFile(ggb_filename,'w') as f:
        f.write("geogebra.xml")
    

if __name__=="__main__":
    
    # Demonstration code
    filename="3vert_on_3face_4verts.ggb"
    G=geogebra_to_graph(filename)
    p=geogebra_graph_plot(G)
    p.save(filename+".pdf")
    
    pos=G.get_pos()
    print("vertices:")
    for v in G.vertices(sort=False):
        info=G.get_vertex(v)
        print(v,info,pos[v])
    print("edges:")
    for e in G.edges(sort=False):
        print(e,G.edge_label(e[0],e[1]))
    
    import sage.graphs.generators.smallgraphs
    G=sage.graphs.generators.smallgraphs.PetersenGraph()
    graph_to_geogebra(G,"petersen.ggb")
    
