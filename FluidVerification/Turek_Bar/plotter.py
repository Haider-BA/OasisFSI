import os
count = 1
get = ["Re", "DOF", "T", "dt", "solver", "theta_scheme", "time", "mesh",]

while os.path.exists("./experiments/cfd3/"+str(count)):
    count += 1

for i in range(1,count):
    info = []
    info.append("Case %d" % i)
    counter = 0
    with open("./experiments/cfd3/"+str(i)+"/report.txt", 'r') as f:
        data = f.readlines()

        for line in data:
            words = line.split()
            for i in words:
                if any(i in s for s in get):
                    #print i, words[-1]
                    info.append(i); info.append(words[-1])

            #print words
    print '[%s]' % ', '.join(map(str, info))
    #print count


"""
with open("./experiments/cfd3/"+str(1)+"/report.txt", 'r') as f:
    data = f.readlines()

    for line in data:
        words = line.split()
        for i in words:
            if any(i in s for s in get):
                print i, words[-1]
        #print words

#print count
"""
