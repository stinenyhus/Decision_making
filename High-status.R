n = 1000
z = 4
n_high = n/z

total_influence = 1/z * n #Influence each node has on its neighbors times number of nodes 
high_influence = n_high * z 
left_influence = total_influence - high_influence
new_influence = left_influence / (n-n_high)

new_influence * z

penalty = 1/z - new_influence

(n - n_high) * new_influence + high_influence

