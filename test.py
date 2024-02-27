a= ['p', 'p', 'q', 'r ', 's', 't', 'y', 'v', 'p', 'p', 'q','p', 'p', 'q','p', 'p', 'q']
testdict = dict()
print(testdict)
aset = set(a)
for itf in aset:
    testdict[itf] = 0
for itf in a:
    testdict[itf] += 1
print(testdict)