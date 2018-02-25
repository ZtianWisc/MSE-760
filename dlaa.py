import csv

data = list(csv.reader(open('20160919_R_TG.csv')))
rows = len(data)
cols = len(data[0])
for i in range(0,rows):
	for j in range(cols):
		data[i][j] = float(data[i][j])
print "m/z   mobility intensity"

for i in range(0,rows-1):
		mzij = data[i+1][0] - data[i][0]
		dtij = data[i+1][1] - data[i][1]
		intensity = data[i+1][2] + data[i][2]
		if (abs(mzij) <= 0.0001 and abs(dtij) >= 1.0 and intensity >= 1000):
			print data[i][0], data[i][1], data[i][2], "|", data[i+1][0], data[i+1][1], data[i+1][2]
		