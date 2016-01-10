#rchen40.py
#Richard Chen

from random import randrange
from math import sqrt

class Mean:
    """Object Mean"""
    
    def __init__(self, mean, close_coords):
        """Object Constructor for Mean. Creates a Mean"""
        self.mean = mean
        self.close_coords = close_coords
    
    def getMean(self):
        """Returns the mean"""
        return self.mean

    def getCloseCoords(self):
        """Returns a list of the closest coordinates"""
        return self.close_coords

    def setMean(self, mean):
        """Sets the mean"""
        self.mean = mean

    def appendCloseCoords(self, coordinates):
        """appends a coordinate to the list of closest coordinates of a mean"""
        self.close_coords.append(coordinates)

    def setCloseCoords(self, close_coords):
        """Sets the closest coordinates"""
        self.close_coords = close_coords

    def __str__(self):
        """Prints out Mean with Closest Coords"""
        return 'Mean: {0}\nClosest Coordinates: {1}'.format(self.mean, self.close_coords)

#Calculates the magnitude of the distance between a coordinate and a mean
def distance(coordinates, mean):
    return sqrt((coordinates[0]-mean[0])**2 + (coordinates[1]-mean[1])**2 + (coordinates[2]-mean[2])**2)

#Calculates the magnitude of the distance between a coordinate and a mean
def distance_18_dimensions(coordinates, mean):
    return sqrt((coordinates[1]-mean[0])**2 + (coordinates[2]-mean[1])**2 + (coordinates[3]-mean[2])**2 + (coordinates[4]-mean[3])**2 + (coordinates[5]-mean[4])**2 + (coordinates[6]-mean[5])**2 + (coordinates[7]-mean[6])**2 + (coordinates[8]-mean[7])**2 + (coordinates[9]-mean[8])**2 + (coordinates[10]-mean[9])**2 + (coordinates[11]-mean[10])**2 + (coordinates[12]-mean[11])**2 + (coordinates[13]-mean[12])**2 + (coordinates[14]-mean[13])**2 + (coordinates[15]-mean[14])**2 + (coordinates[16]-mean[15])**2 + (coordinates[17]-mean[16])**2 + (coordinates[18]-mean[17])**2)
	
#Finds the shortest distance (smallest number) in a list
def find_shortest_distance(distance_list):
    shortest = 999
    
    for distance in distance_list:
        if distance < shortest:
            shortest = distance
            
    return distance_list.index(shortest)
    
##################################################################################################################
#Part 1
##################################################################################################################
    
def kmeans(mean_list, coordinates_list):
    only_means = []
    only_new_means = []
    distance_list = []
    count = 1
    for mean in mean_list:
        only_means.append(mean.getMean())
        

    #This while loop continues until the old mean list is no longer the new mean list.
    while True:
        #This for loop assigns every coordinate to the closest mean.
        for coordinates in coordinates_list:
            #This for loop creates a list of distances of a cooordinate to the different means respectively
            for mean in mean_list:
                distance_list.append(distance(coordinates, mean.getMean()))
            
            shortest_distance_index = find_shortest_distance(distance_list)
            mean_list[shortest_distance_index].appendCloseCoords(coordinates)

            distance_list = []

        #Averages all the coordinate points for each mean
        for mean in mean_list:
            x_tot, y_tot, z_tot = 0, 0, 0
            
            #Total XYZ
            for coordinates in mean.getCloseCoords():
                x_tot = x_tot + coordinates[0]
                y_tot = y_tot + coordinates[1]
                z_tot = z_tot + coordinates[2]

            #This takes the average
            if len(mean.getCloseCoords()) == 0:
                only_new_means.append([x_tot, y_tot, z_tot])
            else:
                only_new_means.append([x_tot/len(mean.getCloseCoords()), y_tot/len(mean.getCloseCoords()), z_tot/len(mean.getCloseCoords())])

       
        #if the new mean is equal to the old mean, break from loop
        if only_means == only_new_means:
            break

        #Updates mean_list and only_means, resets only_new_means
        only_means = []
        for mean in mean_list:
            mean.setMean(only_new_means.pop(0))
            only_means.append(mean.getMean())
            mean.setCloseCoords([])
        only_new_means = []

    return mean_list

def part1():
    print('Part 1\n')
    infilename = input("Enter a file name for K-Means Clustering: ")
    clusters = open(infilename, "r")
    clusters_list = clusters.readlines()
    length = len(clusters_list)

    coordinates_list = []
    #This for loop creates a list of coordinates, with each coordinate a list of three floats: x, y, and z
    for coordinates in clusters_list:
        coordinates = coordinates.split()
        coordinates[0] = float(coordinates[0])
        coordinates[1] = float(coordinates[1])
        coordinates[2] = float(coordinates[2])
        coordinates_list.append(coordinates)
        
    for k in range(2,7):
        index_list = []
        mean_list = []
        #This for loop randomly selects "k" number of coordinates as the mean.
        #The mean is stored in the object "Mean", which not only contains the mean, but the closest coordinates
        #associated with that mean
        for i in range(k):
            rand = randrange(0, len(coordinates_list))
            while (rand in index_list) == True:
                rand = randrange(0, length)
            index_list.append(rand)
            mean_list.append(Mean(coordinates_list[rand],[]))
            
        converged_mean_list = kmeans(mean_list, coordinates_list)

        cluster_index = 1
        
        print('K =', k)
        print('Converged Mean')

        outfile = open("k"+str(k)+".dat", "w")
        for mean in converged_mean_list:
            print(mean.getMean())
            for coordinates in mean.getCloseCoords():
                str_coordinates = ''
                for coordinate in coordinates:
                    str_coordinates = str_coordinates + str(coordinate) + '\t'
                str_coordinates = str_coordinates + str(cluster_index)
                print(str_coordinates, file = outfile)
            cluster_index = cluster_index + 1

        print('\n', file = outfile)
        
        cluster_index = 1
        print()

##################################################################################################################
#Part 2
##################################################################################################################

#KMEANS_18_dimensions
def kmeans_18_dimensions(mean_list, coordinates_list):
    only_means = []
    only_new_means = []
    distance_list = []
    count = 1
    
    for mean in mean_list:
        only_means.append(mean.getMean())

    print("Please be patient and allow this program to execute. With K = 50 and 8269 data points, it might take awhile to converge to the final means.")

    #This while loop continues until the old mean list is no longer the new mean list.
    while True:
        #This for loop assigns every coordinate to the closest mean.
        for coordinates in coordinates_list:
                
            #This for loop creates a list of distances of a cooordinate to the different means respectively
            for mean in mean_list:
                distance_list.append(distance_18_dimensions(coordinates, mean.getMean()))
            
            shortest_distance_index = find_shortest_distance(distance_list)
            mean_list[shortest_distance_index].appendCloseCoords(coordinates)

            distance_list = []
            
        #Averages all the coordinate points for each mean
        for mean in mean_list:
            dim_1, dim_2, dim_3, dim_4, dim_5, dim_6, dim_7, dim_8, dim_9, dim_10, dim_11, dim_12, dim_13, dim_14, dim_15, dim_16, dim_17, dim_18 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            
            #Total Dimensions
            for coordinates in mean.getCloseCoords():
                dim_1 = dim_1 + coordinates[1]
                dim_2 = dim_2 + coordinates[2]
                dim_3 = dim_3 + coordinates[3]
                dim_4 = dim_4 + coordinates[4]
                dim_5 = dim_5 + coordinates[5]
                dim_6 = dim_6 + coordinates[6]
                dim_7 = dim_7 + coordinates[7]
                dim_8 = dim_8 + coordinates[8]
                dim_9 = dim_9 + coordinates[9]
                dim_10 = dim_10 + coordinates[10]
                dim_11 = dim_11 + coordinates[11]
                dim_12 = dim_12 + coordinates[12]
                dim_13 = dim_13 + coordinates[13]
                dim_14 = dim_14 + coordinates[14]
                dim_15 = dim_15 + coordinates[15]
                dim_16 = dim_16 + coordinates[16]
                dim_17 = dim_17 + coordinates[17]
                dim_18 = dim_18 + coordinates[18]

            #New Mean
            if len(mean.getCloseCoords()) == 0:
                only_new_means.append([dim_1, dim_2, dim_3, dim_4, dim_5, dim_6, dim_7, dim_8, dim_9, dim_10, dim_11, dim_12, dim_13, dim_14, dim_15, dim_16, dim_17, dim_18])
            else:
                only_new_means.append([dim_1/len(mean.getCloseCoords()), dim_2/len(mean.getCloseCoords()), dim_3/len(mean.getCloseCoords()), dim_4/len(mean.getCloseCoords()), dim_5/len(mean.getCloseCoords()), dim_6/len(mean.getCloseCoords()), dim_7/len(mean.getCloseCoords()), dim_8/len(mean.getCloseCoords()), dim_9/len(mean.getCloseCoords()), dim_10/len(mean.getCloseCoords()), dim_11/len(mean.getCloseCoords()), dim_12/len(mean.getCloseCoords()), dim_13/len(mean.getCloseCoords()), dim_14/len(mean.getCloseCoords()), dim_15/len(mean.getCloseCoords()), dim_16/len(mean.getCloseCoords()), dim_17/len(mean.getCloseCoords()), dim_18/len(mean.getCloseCoords())])

        #if the new mean is equal to the old mean, break from loop
        if only_means == only_new_means:
            for i in range(len(only_means)):
                print("Mean {0}:".format(i+1), only_means[i])
            break

        #Updates mean_list and only_means, resets only_new_means
        print('Update:', count)
        count = count + 1
        only_means = []
        for mean in mean_list:
            mean.setMean(only_new_means.pop(0))
            only_means.append(mean.getMean())
            mean.setCloseCoords([])
        only_new_means = []

    return mean_list

def part2():
    print('Part 2')
    index_list = []
    mean_list = []
    infilename = input("Enter a file name for the human expression data: ")
    clusters = open(infilename, "r")
    clusters_list = clusters.readlines()
    length = len(clusters_list)

    coordinates_list = []
    clone_coordinates_list = []
    #This for loop creates a list of coordinates, with each coordinate a list of 18 floats and the gene ID
    for coordinates in clusters_list:
        clone_coordinates = coordinates.split()
        coordinates = coordinates.split()
        
        for i in range(1,19):
            coordinates[i] = float(coordinates[i])
            clone_coordinates[i] = float(clone_coordinates[i])
            
        coordinates_list.append(coordinates)
        clone_coordinates.pop(0)
        clone_coordinates_list.append(clone_coordinates)


    #randomly selects 50 coordinate points
    for i in range(50):
            rand = randrange(0, len(clone_coordinates_list))
            while (rand in index_list) == True:
                rand = randrange(0, length)
            index_list.append(rand)
            mean_list.append(Mean(clone_coordinates_list[rand],[]))

    #KMEANS
    converged_mean_list = kmeans_18_dimensions(mean_list, coordinates_list)

    #Formats output
    height = 1000
    tab = open("K50.txt", "w")

    label = ""
    for i in range(50):
        label = label + 'k{0}'.format(i+1) + '\t'
        
    print(label, file = tab)
    for i in range(height):
        str_coordinates = ''
        for mean in converged_mean_list:
            try:
                str_coordinates = str_coordinates + mean.getCloseCoords()[i][0] + '\t'
            except(BaseException):
                str_coordinates = str_coordinates + '\t'
        print(str_coordinates, file = tab)   
    print()
    print("The data has been outputted to K50.txt, with each column as a cluster.")

#main
def main():
    part1()
    part2()
main()
