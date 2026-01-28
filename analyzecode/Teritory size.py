import csv
from sklearn.neighbors import KDTree

#file = open('col3D1.csv')
#csvreader = csv.reader(file)
#header = []
#header = next(csvreader)
#rows = []
#for row in csvreader:
#    rows.append(row)
#rows = rows[2:]
#number_of_fish = (len(row) - 1)/18
#for i in range(0,numberOfFish):
    #print(i)
     #[[[1,2],[2,3]],[],[],[],[],[]]
     #todo create list of x,y points for each fish + reduce number of points based of likelihood
fishPointsList = [[[1,1],[1,2],[1,3],[1,10]]]
avg_fish_distances = []
for fishPoints in fishPointsList:
    tree = KDTree(fishPoints)
    points_avg_distances = []
    for fishPoint in fishPoints:
        dist, ind = tree.query([fishPoint], k=3)
        point_dist = sum(dist[0]) / len(dist[0])
        points_avg_distances += [point_dist]
    avg_fish_distance = sum(points_avg_distances) / len(points_avg_distances)
    avg_fish_distances += [avg_fish_distance]
    print(avg_fish_distance)

