# matplotlib inline
import matplotlib.pyplot as plt
import SALib.analyze.morris

def sum_tuples(l):
    return list(sum(x) for x in zip(*l))

Si = SALib.analyze.morris.analyze(problem, np.array(param_values), 
                                          np.array(sum_tuples(list(prov_roads['max_ini_adap_cost']))),
                                     print_to_console=False, grid_jump=2, num_levels=4)

sensitivity = pd.DataFrame.from_dict(Si)
sensitivity['rel'] = sensitivity['mu']/sensitivity['mu'].sum()*100
sensitivity = sensitivity.groupby('names').sum()
sensitivity = sensitivity.T

stats=sensitivity.loc['rel',np.array(sensitivity.columns)].values

fig, ax = plt.subplots(1, 1,figsize=(6,6),subplot_kw=dict(projection='polar')) #

angles=np.linspace(0, 2*np.pi, len(np.array(sensitivity.columns)), endpoint=False)
# close the plot
stats=np.concatenate((stats,[stats[0]]))
angles=np.concatenate((angles,[angles[0]]))

fig=plt.figure()
ax.plot(angles, stats, 'o-', linewidth=2)
ax.set_ylim([0, 50])   
ax.set_yticklabels([])
ax.fill(angles, stats, alpha=0.25)
ax.set_thetagrids(angles * 180/np.pi, np.array(sensitivity.columns))
ax.tick_params(axis='x',labelsize=10,labelcolor='black',color='black') # pad=12

ax.set_title('Sensitivity Analysis', y=1.1,fontweight='bold')
ax.grid(True)