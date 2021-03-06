#Exercise 7
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn import linear_model, preprocessing

data = pd.read_csv("\\Users\\Henrik\\PycharmProjects\\DataAnalytics-hec-ss19\\Exercise_2\\dat\\data.txt",
                   delim_whitespace = True,
                   names = ["BT", "WBC", "DS", "I"], header = None)

BT, WBC, DS, Inf = data.iloc[:, 0], data.iloc[:, 1], data.iloc[:, 2], data.iloc[:, 3]

BT_h = []
BT_s = []
for i in list(range(len(Inf))):
    if Inf[i] == -1:
        BT_h.append(BT[i])
    else:
        BT_s.append(BT[i])

WBC_h = []
WBC_s = []
for i in list(range(len(Inf))):
    if Inf[i] == -1:
        WBC_h.append(WBC[i])
    else:
        WBC_s.append(WBC[i])

DS_h = []
DS_s = []
for i in list(range(len(Inf))):
    if Inf[i] == -1:
        DS_h.append(DS[i])
    else:
        DS_s.append(DS[i])


s1 = plt.scatter(WBC_h, BT_h, color = "b")
s2 =plt.scatter(WBC_s, BT_s, color = "r")
plt.xlabel("White blood cells")
plt.ylabel("Body temperature(C)")
plt.legend((s1, s2), ("Healthy", "Infected"))
plt.show()

s1 = plt.scatter(DS_h, WBC_h, color = "b")
s2 = plt.scatter(DS_s, WBC_s, color = "r")
plt.xlabel("Daily sleeps(hrs)")
plt.ylabel("White blood cells")
plt.legend((s1, s2), ("Healthy", "Infected"))
plt.show()

s1 = plt.scatter(DS_h, BT_h, color = "b")
s2 = plt.scatter(DS_s, BT_s, color = "r")
plt.xlabel("Daily sleeps(hrs)")
plt.ylabel("Body temperature(C)")
plt.legend((s1, s2), ("Healthy", "Infected"))
plt.show()

# b) Drawing the conditional distribution pairs

h1 = plt.hist(WBC_s, density = True, color = "r")
h2 = plt.hist(WBC_h, density = True, color = "b")
legend = ["Infected", "Healthy"]
plt.ylabel("Probability")
plt.xlabel("White blood cells")
plt.title("Probability distribution of white blood cells")
plt.legend(legend)
plt.show()

h1 = plt.hist(BT_s, density = True, color = "r")
h2 = plt.hist(BT_h, density = True, color = "b")
legend = ["Infected", "Healthy"]
plt.ylabel("Probability")
plt.xlabel("Body temperature(C)")
plt.title("Probability distribution of boody temperature(C)")
plt.legend(legend)
plt.show()

h1 = plt.hist(DS_s, density = True, color = "r")
h2 = plt.hist(DS_h, density = True, color = "b")
legend = ["Infected", "Healthy"]
plt.ylabel("Probability")
plt.xlabel("Daily sleeps(hrs)")
plt.title("Probability distribution of daily sleeps(hrs)")
plt.legend(legend)
plt.show()



# d)

corr = data.corr()
corr.style.background_gradient(cmap='coolwarm')
print(corr)
print("-----------")

# e) Plots for information on single decision algorithms

BT, WBC, DS, Inf = data.iloc[:, 0].values, data.iloc[:, 1].values, data.iloc[:, 2].values, data.iloc[:, 3].values
#Plotting the regression for Infected - White blood cells
Inf = Inf.reshape(Inf.size, 1)
poly = preprocessing.PolynomialFeatures(1)
Infp = poly.fit_transform(Inf)

reg = linear_model.LinearRegression()
reg.fit(Infp, WBC)
model = sm.OLS(WBC, Infp)
y_pred = model.predict(model.fit().params)

plt.scatter(Inf, WBC, color = "black")
plt.plot(Inf, reg.predict(Infp), color = "b", lw = 2)
plt.title("Infected vs White blood cells")
plt.xlabel("Infected")
plt.ylabel("White blood cells")
plt.show()

#Plotting the regression for Infected - body tempreature
reg = linear_model.LinearRegression()
reg.fit(Infp, BT)
model = sm.OLS(BT, Infp)
y_pred = model.predict(model.fit().params)

plt.scatter(Inf, BT, color = "black")
plt.plot(Inf, reg.predict(Infp), color = "b", lw = 2)
plt.title("Infected vs body temperature(C)")
plt.xlabel("Infected")
plt.ylabel("body temperature(C)")
plt.show()

#Plotting the regression for Infected - daily hours of sleep

reg = linear_model.LinearRegression()
reg.fit(Infp, DS)
model = sm.OLS(DS, Infp)
y_pred = model.predict(model.fit().params)

plt.scatter(Inf, DS, color = "black")
plt.plot(Inf, reg.predict(Infp), color = "b", lw = 2)
plt.title("Infected vs Daily sleep(hrs)")
plt.xlabel("Infected")
plt.ylabel("Daily hours of sleep")
plt.ylim(2,10)
plt.show()

