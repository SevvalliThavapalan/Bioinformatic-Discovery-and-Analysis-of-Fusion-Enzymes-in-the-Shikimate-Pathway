import matplotlib.pyplot as plt
import numpy as np

enzyme = [49.83, 42.127, 4.118, 1.783, 0.769, 0.556, 0.491, 0.327]
pos = ["Others", "AroB", "p dehydra", "AroQ", "AroD", "AroA", "AroC", "AroG" ]
fig = plt.figure(figsize=(15,40))
ax = fig.add_subplot(111)
ax.bar(x=pos, height=enzyme, color='#1f77b4')
plt.xticks(rotation=45, fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel("Percent", fontsize=20)
plt.title("AroK - Right enzymes",fontsize=20 )


 # 50
plt.show()
