
# loss Function Next:
loss = 0
N = len(X)
for k in range(0,N):
    loss = loss + (grid_Yhat[k] - Y[k]) * (grid_Yhat[k] - Y[k])
lossupdated = math.sqrt(loss/N)
print lossupdated

#Check for overfitting
loss = 0
N = len(X)
for k in range(N * 4/5,N):
    loss = loss + (grid_Yhatvalidation[k] - Y[k]) * (grid_Yhatvalidation[k] - Y[k])
lossupdated = math.sqrt(loss/N)
print lossupdated