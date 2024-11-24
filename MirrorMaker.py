import numpy as np
import csv
from scipy import interpolate
import warnings

## Defining functions.

def read_vessel_file(filename, delimiter=" "):
    """
    This function reads the vessel file and transforms it into a form which can be used to make the mirror files.
    
    filename: the filename of the vessel file.
    delimiter: the delimiter used in the vessel file. By default a space.
    """
    
    ## Import the data
    data = []
    with open(filename, 'r') as file:
        # Create a CSV reader object, specifying the delimiter as a tab
        reader = csv.reader(file, delimiter=delimiter)

        # Read the data row by row
        for row_reader in reader:
            row = []
            for item in row_reader:
                if item!="":# Due to how the reader works a lot of items in the row list will be empty, ignore those.
                    try:
                        row.append(float(item))
                    except ValueError:# Only accept numbers, if it is a string, ignore it.
                        pass
            
            data.append(row)
    
    ## Transform the data into a format which can be better used.
    transformed_data = []
    transformed_row=[data[2][0]]
    
    for i in range(3,len(data)):# We skip the first line because that is already done while generating the transformed data and row list.
        if len(data[i])==1:# Make a row/list for each phi, each row contains sub-list with the r and z coordinates.
            transformed_data.append(transformed_row)
            transformed_row=[data[i][0]]
        else:
            transformed_row.append(data[i])
        
    return transformed_data


def make_mirrors(vessel_data, mirror_locations):
    """
    This function determines the coordinates of the mirror based on the vessel data and specified locations. The specifies phi range is rounded to the closest matching phi in the vessel file. This function does't work near the top and bottom of the reactor. Only near the sides.
    
    vessel_data: is the vessel data (as exported by read vesel file. Each (phi) segment is a list element with the r and z coordinates as a sublist for each phi.
    mirror_locations: is a list consisting of dictionairies. Each dictionairy specifies if it is inboard or outboard, the phi and z location (from low to high) and can optionally (but highly recommended) specify margin (same units as vessel file) and the number of subdivsions in z. Example: [{"location":"inboard", "phi range":[0,4], "z range":[-45, 20], "number of z points":10, "margin":1}, {"location":"outboard", "phi range":[0,4], "z range":[-45, -20]}]
    """
    mirrors = []# Array for storing each mirror in.
    phis = np.array([segment[0] for segment in (vessel_data)])# Get all the phi values present in the vessel file.
    
    for num in range (len(mirror_locations)):# Loop through each mirror location
        ## Data extraction (if vallues are missing set default values
        phi_range = mirror_locations[num]["phi range"]
        z_range = mirror_locations[num]["z range"]
        try:
            num_z_points = mirror_locations[num]["number of z points"]
        except KeyError:
            num_z_points = 2
            warnings.warn("Number of z points of mirror " + str(num) + " not specified. Default value of 2 used.")
        
        try:
            margin = mirror_locations[num]["margin"]
        except KeyError:
            margin = 0
            warnings.warn("Margin of mirror " + str(num) + " not specified. Default value of 0 used, which can cause problems.")
        
        if mirror_locations[num]["location"].lower()=="inboard":
            location=1
        elif mirror_locations[num]["location"].lower()=="outboard":
            location=-1
        else:
            raise ValueError("Location doesn't specify if the mirrors should be 'Inboard' or 'Outboard'.")
        
        index_phi_start = (np.abs(phis - phi_range[0])).argmin()# Determine the index of the relevant phies
        index_phi_end = (np.abs(phis - phi_range[1])).argmin()
        
        if index_phi_end==index_phi_start:
            index_phi_end = index_phi_end+1
            warning.warn("Phi range is smaller than phi resolution or is outside of the phi range of the vessel file. The cod automatically extended the phi range of the mirrror to prevent errors.")
        
        
        polygons=[]# List to store each polygon of the mirror in. Each polygon is a rectangle.
        for segment_num in range(index_phi_start, index_phi_end):
            segment = vessel_data[segment_num].copy()
            phi = segment.pop(0)
            major_radius = np.mean([point[0] for point in segment])# Calculate the major radius to use to detemine what is the inboard and outboard side.
            
            ## Get the inboard/outboard coordinates of this segment and the next segment (both are needed to make the rectangle).
            r_coordinates = []
            z_coordinates = []
            for i in range(len(segment)):
                if segment[i][0]**location<major_radius**location:# This is to determine the inboard or outboard side. Uses the property that if x<y then 1/x>1/y.
                    r_coordinates.append(segment[i][0])
                    z_coordinates.append(segment[i][1])
            
            interp = interpolate.interp1d(z_coordinates, r_coordinates, fill_value='extrapolate')
                
            next_segment = vessel_data[segment_num+1].copy()
            next_phi = next_segment.pop(0)
            next_major_radius = np.mean([point[0] for point in segment])
            
            next_r_coordinates = []
            next_z_coordinates = []
            for i in range(len(next_segment)):
                if next_segment[i][0]**location<next_major_radius**location:# This is to determine the inboard or outboard side. Uses the property that if x<y then 1/x>1/y.
                    next_r_coordinates.append(next_segment[i][0])
                    next_z_coordinates.append(next_segment[i][1])
            
            next_interp = interpolate.interp1d(next_z_coordinates, next_r_coordinates, fill_value='extrapolate')
            
            ## Now use interpolation to determine the r at the specified z values. Then do the same for the next segmant. Each z point uses a seperate polygon/rectangle.
            for z_num in range(num_z_points):
                polygon=[]
                
                delta = (z_range[1] - z_range[0])/num_z_points
                z_polygon = np.array([z_range[0]+z_num*delta, z_range[0]+(z_num+1)*delta])#np.linspace(z_range[0], z_range[1], num=num_z_segments)
                r_polygon = interp(z_polygon)
                phi_polygon = phi*np.ones(2)
                
                next_z_polygon = np.array([z_range[0]+(z_num+1)*delta, z_range[0]+z_num*delta])#np.linspace(z_range[1], z_range[0], num=num_z_segments)
                next_r_polygon = next_interp(next_z_polygon)
                next_phi_polygon = next_phi*np.ones(2)
                
                z_polygon = np.concatenate((z_polygon, next_z_polygon))
                r_polygon = np.concatenate((r_polygon, next_r_polygon))
                phi_polygon = np.concatenate((phi_polygon, next_phi_polygon))
                
                for i in range(len(z_polygon)):# Add the margin and then calculate the x and y coordinates and save it as a polygon vertex
                    radius_total_initial = np.sqrt(z_polygon[i]**2+r_polygon[i]**2)
                    r = r_polygon[i] + location*margin#*r_polygon[i]/radius_total_initial
                    z = z_polygon[i]# - location*margin*z_polygon[i]/radius_total_initial
                    #radius_total = np.sqrt(r**2+z**2)
                    
                    x = r*np.cos(phi_polygon[i]*np.pi/180)
                    y = r*np.sin(phi_polygon[i]*np.pi/180)
                    vertex = [x, y, z]
                    
                    polygon.append(vertex)
                
                polygons.append(polygon)
            
        mirrors.append(polygons)
        
    return mirrors

def make_mirror_file(mirrors, filename, conversion_to_mm=0.1, comment="", tolerance=0.005):
    """
    This function converts the mirros provided by the make_mirrors into a file.
    
    mirrors: a list of depth 4 which describes the mirror(s) as generated by the make_mirrors function.
    conversion_to_mm: is the conversion factor to the unit of the mirror file: milimeters. By default is 0.1 (centimets per milimeter).
    comment: is a comment which is added to the XML file. It is recommended to use this to document which vesselfile was used and what settings for the mirror.
    """
    file = open(filename, 'w')
    
    file.write("<?xml version=\"1.0\" ?>\n")
    file.write("<!--"+comment+"-->\n")
    file.write("<mirror id=\"CustomMirror\" description=\"ECRH mirror\">\n")
    
    for surface_num in range(len(mirrors)):
        surface = mirrors[surface_num]
        file.write("\t<surface id = \""+str(surface_num)+"\" name = \""+str(surface_num)+"\">\n")
        
        for polygon in surface:
            file.write("\t\t<polygon unit = \"mm\" >\n")
            for vertex in polygon:
                file.write("\t\t\t<vertex tolerance = \""+ str(tolerance)+"\">"+str(vertex[0]/conversion_to_mm) + " " + str(vertex[1]/conversion_to_mm) + " " + str(vertex[2]/conversion_to_mm) + "</vertex>\n")
            
            file.write("\t\t</polygon>\n")
        
        file.write("\t</surface>\n")
    
    file.write("</mirror>\n")
    file.close() 

## Running the script
filepath = 'wout_squid_20230921_v1_vessel_hightotal_res.txt'
mirror_locations = [{"location":"inboard", "phi range":[85,95], "z range":[-10, 10], "number of z points":10, "margin":1}, {"location":"outboard", "phi range":[175,185], "z range":[-10, 10], "number of z points":10, "margin":1}]
#, {"location":"outboard", "phi range":[72-2,72+4], "z range":[-20, -11], "number of z points":10, "margin":1}, {"location":"outboard", "phi range":[72+68,72+74], "z range":[-10, 10], "number of z points":10, "margin":1}]

vessel_data = read_vessel_file(filepath)
mirrors = make_mirrors(vessel_data, mirror_locations)
make_mirror_file(mirrors, "Squid_mirror_own.xml", comment=filepath+str(mirror_locations))