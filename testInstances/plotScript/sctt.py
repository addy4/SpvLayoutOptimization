import matplotlib.pyplot as plt 
import subprocess as sp 
  
# dataset-1 
x1 = [89, 43, 36, 36, 95, 10,  
      66, 34, 38, 20] 
  
y1 = [21, 46, 3, 35, 67, 95,  
      53, 72, 58, 10] 
  
# dataset2 
x2 = [26, 29, 48, 64, 6, 5, 
      36, 66, 72, 40] 
  
y2 = [26, 34, 90, 33, 38,  
      20, 56, 2, 47, 15]

x11 = [] 
x22 = [] 

file_name_med = "median_points4" 
file_name_panels = "panels_points_11"  

med_lines = sp.getoutput("wc -l median_points4 | awk '{print $1}' ") 
panels_lines = sp.getoutput("wc -l panels_points_11 | awk '{print $1}' ") 

panels_lines = int(panels_lines) 
med_lines = int(med_lines) 

for i in range(panels_lines) : 
	command = "head -{0} {1} | tail -1".format(i+1, "panels_points_11")    
	print(command) 
	command_for_N = command + " | awk "  +  "'{" + " print $1 " + "}'"    
	N = sp.getoutput(command_for_N) 
	N = int(N) 
	x11.append(N)
	command_for_M = command + " | awk "  +  "'{" + " print $2 " + "}'"
	M = sp.getoutput(command_for_M) 
	M = int(M) 
	x22.append(M) 

y11 = [] 
y22 = [] 

for j in range(med_lines) :
        command = "head -{0} {1} | tail -1".format(j+1, "median_points4")  
        print(command)
        command_for_N = command + " | awk "  +  "'{" + " print $3 " + "}'"
        N = sp.getoutput(command_for_N)
        N = int(N)
        y11.append(N)
        command_for_M = command + " | awk "  +  "'{" + " print $6 " + "}'"
        M = sp.getoutput(command_for_M)
        M = int(M)
        y22.append(M) 
	

print(x11) 
print(x22) 


plt.scatter(x11, x22, c ="pink",  
            linewidths = 2,  
            marker ="s",  
            edgecolor ="green") 
  
plt.scatter(y11, y22, c ="yellow", 
            linewidths = 2, 
            marker ="^",  
            edgecolor ="red",s=75) 
  
plt.xlabel("X-axis") 
plt.ylabel("Y-axis") 
plt.show()
#plt.savefig("Sparse_images/fig_panelPoints10") 
