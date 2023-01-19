import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme()
beta = 1.66        #declare given constant variables
R_not = 3.65
gamma = 1.66/R_not
N = 763

r,I = [0],[1/N] #store S,I and R values
s =[1-I[0]-r[0]]

#defining functions for S,I,R respectively
def dsdt(s, I):
  return -beta*s*I
def dIdt(s,I):
  return beta*s*I - gamma*I
def drdt(I):
  return gamma*I

#defining the algorithm for SIR model using 4th order Runge-Kutta Method
def RK4SIR(s,I,r, h = 0.25,time = [0]):
  i=0
  #find values for S[i]s and I[i]s
  while time[i] < 30: #calculating for 30 days
    sk1 = dsdt(s[i],I[i])
    Ik1 = dIdt(s[i],I[i])

    sk2 = dsdt(s[i] + h*sk1/2, I[i] + h*Ik1/2)
    Ik2 = dIdt(s[i] + h*sk1/2, I[i] + h*Ik1/2)

    sk3 = dsdt(s[i] + h*sk2/2, I[i] + h*Ik2/2)
    Ik3 = dIdt(s[i] + h*sk2/2, I[i] + h*Ik2/2)

    sk4 = dsdt(s[i] + h*sk3/2, I[i] + h*Ik3/2)
    Ik4 = dIdt(s[i] + h*sk3/2, I[i] + h*Ik3/2)

    si = s[i]+h*(sk1 + 2*sk2 + 2*sk3 + sk4)/6
    Ii = I[i]+h*(Ik1 + 2*Ik2 + 2*Ik3 + Ik4)/6
    s.append(si)
    I.append(Ii)
    time.append(time[i] + h)
    i+=1
  
  #Since S+I+R = 1 => R = 1-S-I
  for i in range(1,len(s)):
    r.append(1-s[i]-I[i])
  return s,I,r,time

#Function definition to plot the results
def plot(s,I,r,time):
  fig, ax = plt.subplots()
  ax = sns.lineplot(x=time, y=s, label = "S")

  ax1 = sns.lineplot(x=time, y=I, label = "I")

  ax3 = sns.lineplot(x = time, y = r, label = "R")

  plt.xlabel("Number of Days")
  plt.title("Spread of Influenza Virus")
  plt.show()

s,I,r,time = RK4SIR(s,I,r) #returns result in respective variables


print("Maximum Value of I: ", max(I)) #returns max value = 0.375117
TIME = time[I.index(max(I))] #time(days) to reach max value of I
print("Time to reach maximum value:", TIME, "Days (", TIME*24, "hours )")
print()

#interpolate by observing both array I and n values to find 1% infection Time(in hours) to reach 1% infection.
x = (((time[76] - time[77])/(I[76]-I[77]))*(0.01*max(I) - I[76])) + time[76] 
print("Time to reach 1% infection:","%.2f" % x , "Days (", "%.2f" % (x*24), "hours )")


x = x-TIME
print("Time taken after peak to reach 1% infection:","%.2f" % x , "Days (", "%.2f" % (x*24), "hours )")


plot(s,I,r,time) #plot the results