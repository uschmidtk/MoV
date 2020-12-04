'''
    File: experiments.py
    Author: Markus Brill
            Ulrike Schmidt-Kraepelin (u.schmidt-kraepelin@tu-berlin.de)
            Warut Suksompong
    Date:   Dec. 4, 2020 
'''
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import itertools
import random
import math

import generate_profiles as pl

####################################################
################# CP FUNCTIONS #####################
####################################################

def cp_scores(T):
    scores = [x[1] for x in sorted(T.out_degree(T.nodes()))]
    return scores
    
def cp_winners(T):
    max_degree = max(cp_scores(T))
    copeland_winners = [x[0] for x in list(T.out_degree(T.nodes())) if x[1]==max_degree]
    return copeland_winners
    
def second_highest(list):
    if type(list)==list:
        max_element = max(list)
        list = [x for x in list if not x == max_element]
        return max(list)
    else:
        return list[0]
    
    
def cp_runner_ups(T):
    if not set(cp_winners(T)) == set(T.nodes()):
        runner_degree = second_highest(cp_scores(T))
        runner_ups = [x[0] for x in list(T.out_degree(T.nodes())) if x[1]==runner_degree]
        return runner_ups
    else:
        return []
        
def cp_mov(T,i,winners,runner_ups,scores):
    if i in winners:
        if len(winners) > 1:
            return 1
        else:
            if set(runner_ups).intersection(set(T.neighbors(i))) != set():
                return max([max(scores)-second_highest(scores),1])
            else:
                return max(scores)-second_highest(scores)+1
    else:
        if set(winners).intersection(set(T.predecessors(i)))==set(winners) and len(winners) < max(scores) - T.out_degree(i):
            return -(max(scores) - T.out_degree(i))+1
        else:
            return -(max(scores) - T.out_degree(i))



####################################################
################# UC FUNCTIONS #####################
####################################################

def double_loop(S):
    for i in range(1,len(S)+1):
        for add in itertools.combinations(S,i):
            yield i,add
            
def in_uc(T,s):
    zero_step = set([s])
    one_step = set(T.neighbors(s))
    two_step = set(nx.algorithms.boundary.node_boundary(T,T.neighbors(s)))
    if zero_step.union(one_step).union(two_step) == set(T.nodes()): return True
    else: return False
            
def uc_cut(T,s,t):
    S = T.copy()
    S.remove_nodes_from(set(T.nodes()) - set(list(T.neighbors(s))) - set([s,t]))
    cut2 = nx.minimum_cut(S,s,t)
    cut_value2 = cut2[0]
    return cut_value2
        
def uc_mov(T,i):
    if(in_uc(T,i)):
        min_uc_value = len(T.edges())
        other_nodes = [x for x in T.nodes() if x != i]
        for j in other_nodes:
            uc_value = uc_cut(T,i,j)
            if uc_value <= min_uc_value:
                min_uc_value = uc_value
        return min_uc_value
    else:
        zero_step = set([i])
        neighbors = set(T.neighbors(i))
        non_neighbors = set(T.nodes) - set([i]) - set(T.neighbors(i))
        for j, add in double_loop(non_neighbors):
            one_step = neighbors.union(add)
            two_step = set(nx.algorithms.boundary.node_boundary(T,one_step))
            if zero_step.union(one_step).union(two_step) == set(T.nodes()):
                return -j
        
    
    
####################################################
################# 3KINGS FUNCTIONS #################
####################################################
    
    
def in_3kings(T,s):
    zero_step = set([s])
    one_step = set(T.neighbors(s))
    two_step = set(nx.algorithms.boundary.node_boundary(T,one_step))
    three_step = set(nx.algorithms.boundary.node_boundary(T,two_step))
    if zero_step.union(one_step).union(two_step).union(three_step) == set(T.nodes()): return True
    else: return False
    
def king_cut(T,u,v):
    U=set(T.neighbors(u)) - set([v])
    V=set(T.predecessors(v)) - set([u])
    
    G=nx.DiGraph()
    G.add_nodes_from([u])
    G.add_nodes_from([v])
    G.add_nodes_from(U)
    G.add_nodes_from([x+n for x in V])
    G.add_edges_from(nx.algorithms.edge_boundary(T,[u],U))
    G.add_edges_from([(x,x+n) for x in U.intersection(V)])
    G.add_edges_from([(x+n,v) for x in V])
    G.add_edges_from([(x,y+n) for (x,y) in nx.algorithms.edge_boundary(T,U,V)])

    for s,t in G.edges():
        G[s][t]['capacity'] = 1
        
    cut = nx.minimum_cut(G,u,v)
    cut_value = cut[0]

    if (u,v) in list(T.edges()):
        cut_value = cut_value + 1
    return cut_value
    
def king_mov(T,i):
    if in_3kings(T,i):
        min_king_value = len(T.edges())
        other_nodes = [x for x in T.nodes() if x != i]
        for j in other_nodes:
            king_value = king_cut(T,i,j)
            if king_value <= min_king_value:
                min_king_value = king_value
        return min_king_value
    else:
        return -1
    

    
    
####################################################
################# TC FUNCTIONS #####################
####################################################

def in_tc(T,s):
    for t in T.nodes():
        if not nx.algorithms.tournament.is_reachable(T,s,t):
            return False
    return True

def tc_cut(T,s,t):
    cut = nx.minimum_cut(T,s,t)
    cut_value = cut[0]
    return cut_value
    
    
def tc_mov(T,i):
    if(in_tc(T,i)):
        min_tc_value = len(T.edges())
        other_nodes = [x for x in T.nodes() if x != i]
        for j in other_nodes:
            tc_value = tc_cut(T,i,j)
            if tc_value <= min_tc_value:
                min_tc_value = tc_value
        return min_tc_value
    else:
        return -1

####################################################
################# CALCULATE MOVS ###################
####################################################

def movs(T,cp=True,uc=True,kings=True,tc=True,cp_score=True):
    movs = pd.DataFrame()
    
    if uc or tc:
        for a,b in T.edges():
            T[a][b]['capacity'] = 1
            
    if cp:
        cp_winner = cp_winners(T)
        cp_score = cp_scores(T)
        runner_up = cp_runner_ups(T)
            
    
    for v in T.nodes():
        if cp:
            movs.at[v,'cp_mov'] = cp_mov(T,v,cp_winner,runner_up,cp_score)
            movs.at[v,'cp'] = int(v in cp_winners(T))
        if uc:
            movs.at[v,'uc_mov'] = uc_mov(T,v)
            movs.at[v,'uc'] = int(in_uc(T,v))
        if kings:
            movs.at[v,'kings_mov'] = king_mov(T,v)
            movs.at[v,'kings'] = int(in_3kings(T,v))
        if tc:
            movs.at[v,'tc_mov'] = tc_mov(T,v)
            movs.at[v,'tc'] = int(in_tc(T,v))
        if cp_score:
            movs.at[v,'cp_score'] = T.out_degree(v)
    return movs

####################################################
################# CREATE TOURNAMENTS ###############
####################################################
        
def condorcet_tournament(n,m,p):
    all_edges=[]
    
    for i in range(n):
        for s in itertools.combinations(range(m),2):
            s=list(s)
            if s[0]!=s[1]:
                coin = np.random.rand()
                if ((s[0] < s[1]) and (coin <= p)) or ((s[1] < s[0]) and (coin > p)):
                    all_edges.append((s[0],s[1]))
                else:
                    all_edges.append((s[1],s[0]))
                    
    edge_count = pd.Series(all_edges).value_counts()
    agg_edges = list(edge_count[edge_count >int(n/2)].index)

    T = nx.DiGraph()
    T.add_nodes_from(range(m))
    T.add_edges_from(agg_edges)
    if nx.tournament.is_tournament(T):
        return T
    else: print("There has been a mistake, this is not a tournament!")


def condorcet_tournament_direct(m,p):
    all_edges=[]
    for s in itertools.combinations(range(m),2):
        s = list(s)
        if s[0]!=s[1]:
            coin = np.random.rand()
            if ((s[0] < s[1]) and (coin <= p)) or ((s[1] < s[0]) and (coin > p)):
                all_edges.append((s[0],s[1]))
            else:
                all_edges.append((s[1],s[0]))
    T = nx.DiGraph()
    T.add_nodes_from(range(m))
    T.add_edges_from(all_edges)
    if nx.tournament.is_tournament(T):
        return T
    else: print("There has been a mistake, this is not a tournament!")

def impartial_culture(n,m):
    all_edges=[]
    for i in range(n):
        order = list(np.random.permutation(range(m)))
        for s in itertools.combinations(range(m),2):
            s=list(s)
            if s[0]!=s[1]:
                if (order.index(s[0]) < order.index(s[1])):
                    all_edges.append((s[0],s[1]))
                else:
                    all_edges.append((s[1],s[0]))
                    
    edge_count = pd.Series(all_edges).value_counts()
    agg_edges = list(edge_count[edge_count>int(n/2)].index)

    T = nx.DiGraph()
    T.add_nodes_from(range(m))
    T.add_edges_from(agg_edges)
    if nx.tournament.is_tournament(T):
        return T
    else: print("There has been a mistake, this is not a tournament!")
    
    
def mallows(n,m,phi):
    candmap = {i:i for i in range(m)}
    rankmapcounts = pl.gen_mallows(n,candmap,[1],[phi],[list(range(m))])
    all_edges = []
    for i in range(len(rankmapcounts[1])):
        for s in itertools.combinations(range(m),2):
            if s[0]!=s[1]:
                if rankmapcounts[0][i][s[0]] < rankmapcounts[0][i][s[1]]:
                    for j in range(rankmapcounts[1][i]):
                        all_edges.append((s[0],s[1]))
                else:
                    for j in range(rankmapcounts[1][i]):
                        all_edges.append((s[1],s[0]))
                        
    edge_count = pd.Series(all_edges).value_counts()
    agg_edges = list(edge_count[edge_count>int(n/2)].index)

    T = nx.DiGraph()
    T.add_nodes_from(range(m))
    T.add_edges_from(agg_edges)
    if nx.tournament.is_tournament(T):
        return T
    else: print("There has been a mistake, this is not a tournament!")
    
def urn(n,m,replace):
    candmap = {i:i for i in range(m)}
    rankmapcounts = pl.gen_urn_strict(n,replace,candmap)
    #print(rankmapcounts)
    all_edges = []
    for i in range(len(rankmapcounts[1])):
        for s in itertools.combinations(range(m),2):
            if s[0]!=s[1]:
                if rankmapcounts[0][i][s[0]] < rankmapcounts[0][i][s[1]]:
                    for j in range(rankmapcounts[1][i]):
                        all_edges.append((s[0],s[1]))
                else:
                    for j in range(rankmapcounts[1][i]):
                        all_edges.append((s[1],s[0]))
                        
    edge_count = pd.Series(all_edges).value_counts()
    agg_edges = list(edge_count[edge_count>int(n/2)].index)
    #print(edge_count)

    T = nx.DiGraph()
    T.add_nodes_from(range(m))
    T.add_edges_from(agg_edges)
    if nx.tournament.is_tournament(T):
        return T
    else: print("There has been a mistake, this is not a tournament!")
    
    
####################################################
# Test Conj. 'Non-Winner UC-MOV' takes few values ##
# Just for documentation, these experiments were  ##
# only mentioned as a side remark and not in      ##
# detail                                          ##
####################################################


if False:
    n_list = [5,10,15,20,25,30]
    sample_size = 100
    model_list = ['random','condorcet_direct','impartial','mallows','urn']



    for n in n_list:
        non_winners = pd.DataFrame(columns=['i','model'])
        non_winners = non_winners.set_index(['i','model'])
        
        for model in model_list:
            for i in range(sample_size):
                if model == 'random':
                    T=nx.algorithms.tournament.random_tournament(n)
                if model == 'condorcet_voters':
                    T=condorcet_tournament(51,n,0.55)
                if model == 'condorcet_direct':
                    T=condorcet_tournament_direct(n,0.55)
                if model == 'impartial':
                    T=impartial_culture(51,n)
                if model == 'mallows':
                    T=mallows(51,n,0.95)
                if model == 'urn':
                    T=urn(51,n,10)
                df = movs(T,cp=False,uc=True,kings=False,tc=False,cp_score=False)
                n_loosers = n - df['uc'].sum()
                if(n_loosers>0):
                    fractions = df[df['uc_mov']<0]['uc_mov'].value_counts()/n_loosers
                    fractions.index = fractions.index*(-1)
                    ub = int(math.log2(n)) + 1
                    add_to = pd.Series(index=list(range(1,ub+1)),data=np.zeros(ub))
                    fractions = (fractions + add_to).fillna(0)
                    fractions = fractions.cumsum()
                    for ind in list(fractions.index):
                        non_winners.at[(i,model),ind] = fractions[ind]
        
        agg_non = non_winners.groupby('model').mean()
        
        plt.figure()
        agg_non.plot(kind='bar')
        plt.savefig(str(n) + '-non-winners-uc')
        plt.close()
        
    
    
    
if False:
                    
    n_list = [5,10,15,20,25,30]
    sample_size = 100
    model_list = ['random','condorcet_direct','impartial','mallows','urn']


    non_winners = pd.DataFrame(columns=['i','n'])
    non_winners = non_winners.set_index(['i','n'])

    for n in n_list:
        for model in model_list:
            for i in range(sample_size):
                if model == 'random':
                    T=nx.algorithms.tournament.random_tournament(n)
                if model == 'condorcet_voters':
                    T=condorcet_tournament(51,n,0.55)
                if model == 'condorcet_direct':
                    T=condorcet_tournament_direct(n,0.55)
                if model == 'impartial':
                    T=impartial_culture(51,n)
                if model == 'mallows':
                    T=mallows(51,n,0.95)
                if model == 'urn':
                    T=urn(51,n,10)
                df = movs(T,cp=False,uc=True,kings=False,tc=False,cp_score=False)
                n_loosers = n - df['uc'].sum()
                if(n_loosers>0):
                    unique = len(df[df['uc_mov']<0]['uc_mov'].value_counts().index)
                    non_winners.at[(i,n),model] = unique
        
    print(non_winners)
    agg_non = non_winners.groupby('n').mean()
    print(agg_non)

    plt.figure()
    agg_non.plot(kind='bar')
    plt.savefig('non-winners-uc')
    plt.close()
            


####################################################
################# START EXPERIMENTS ################
####################################################

sample_size = 100 #original sample size
#sample_size = 10 #recommended list for first test run
n_list = [5,10,15,20,25,30] #list from original experiments
#n_list = [3,5,6] #recommended list for first test run
model_list = ['random','condorcet_voters','condorcet_direct','impartial','mallows','urn'] #full list of used models
#model_list = ['random'] #recommended for first test run
run = True #decides whether experiments run
zoomin = False #if set to true produces seperate plots for each size of n,target value and model


title = {'random':'Uniform Random', 'condorcet_voters':'Condorcet Noise via Voters p=0.55, n=51','condorcet_direct':'Condorcet Noise p=0.55','impartial':'Impartial Culture', 'mallows': 'Mallows (phi = 0.95)','urn':'Urn (alpha=10)'}

if run:

    for model in model_list:

        maax = pd.DataFrame(columns=['i','n'])
        maax = maax.set_index(['i', 'n'])
        unique = pd.DataFrame(columns=['i','n'])
        unique = unique.set_index(['i', 'n'])

        for n in n_list:
            for i in range(sample_size):
                if model == 'random':
                    T=nx.algorithms.tournament.random_tournament(n)
                if model == 'condorcet_voters':
                    T=condorcet_tournament(51,n,0.55)
                if model == 'condorcet_direct':
                    T=condorcet_tournament_direct(n,0.55)
                if model == 'impartial':
                    T=impartial_culture(51,n)
                if model == 'mallows':
                    T=mallows(51,n,0.95)
                if model == 'urn':
                    T=urn(51,n,10)
                df = movs(T)
                for col in list(df.columns):
                    maax.at[(i,n),col] = len(df[col][df[col] == df[col].max()])
                for col in ['cp_mov', 'cp_score', 'uc_mov', 'kings_mov', 'tc_mov']:
                    unique.at[(i,n),col] = len(df[col].unique())
                    
            if zoomin:
                for col in list(maax.columns):
                    All = slice(None)
                    maaxi = maax.loc[(All,n),col]
                    plt.figure()
                    maaxi.plot(kind='hist',bins=n, title= col + ' - Histogram of Maximum Equivalence Class - ' + title[model],density=1)
                    plt.savefig(str(n) + '-' + col + '-hist-max-' + model)
                    plt.close()
                    
                for col in ['cp_mov', 'cp_score', 'uc_mov', 'kings_mov', 'tc_mov']:
                    unique2 = unique.loc[(All,n),col]
                    plt.figure()
                    unique2.plot(kind='hist',bins=n, title= col + ' - Histogram of Unique Values - ' + title[model],density=1)
                    plt.savefig(str(n) + '-' + col + '-hist-unique-' + model)
                    plt.close()
                    
                

        maax.to_csv('max-' + model + '.csv')
        agg_maax = maax.groupby('n').mean()
        agg_maax.to_csv('agg_max-' + model + '.csv')
        
        ind_str = np.array(['n = ' + str(x) for x in n_list])
        ind = np.arange(len(n_list))
        width = 7/100
        fig, ax = plt.subplots()
        hfont = {'fontname':'Times'}
        
        rects1 = ax.bar(ind - 0.33, np.array(agg_maax['cp_mov']), width, label='CO-mov',color='#2961A6')
        rects2 = ax.bar(ind - 0.33 + width, np.array(agg_maax['cp']), width, label='CO',color='#76A8D7')
        rects3 = ax.bar(ind - 0.33 + 8/3*width, np.array(agg_maax['uc_mov']), width, label='UC-mov',color='#B02417')
        rects4 = ax.bar(ind - 0.33 + 11/3*width, np.array(agg_maax['uc']), width, label='UC',color='#F19393')
        rects5 = ax.bar(ind - 0.33 + 16/3*width, np.array(agg_maax['kings_mov']), width, label='3-kings-mov',color='#058E41')
        rects6 = ax.bar(ind - 0.33 + 19/3*width, np.array(agg_maax['kings']), width, label='3-kings',color='#60C87B')
        rects7 = ax.bar(ind - 0.33 + 24/3*width, np.array(agg_maax['tc_mov']), width, label='TC-mov',color='#E5AB14')
        rects8 = ax.bar(ind - 0.33 + 27/3*width, np.array(agg_maax['tc']), width, label='TC',color='#F7E091')
        
        # Add some text for labels, title and custom x-axis tick labels, etc.
        #ax.set_ylabel('Average Size')
        ax.set_title(title[model],loc='right',size=14,fontname='Times')
        ax.set(ylim=[0,max(n_list)+1])
        ax.set_xticks(ind)
        ax.set_xticklabels(ind_str)
        ax.tick_params(labelsize=13)
        ax.minorticks_on()
        ax.legend(loc='upper left')
        plt.savefig('max-' + model + '2',dpi=500)
        plt.close()
        
        unique.to_csv('unique-' + model + '.csv')
        agg_unique = unique.groupby('n').mean()
        agg_unique.to_csv('agg_unique-' + model + '.csv')
        
        width = 0.14
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind - 0.27, np.array(agg_unique['cp_mov']), width, label='CO-mov',color='#2961A6')
        rects3 = ax.bar(ind - 0.27 + width, np.array(agg_unique['uc_mov']), width, label='UC-mov',color='#B02417')
        rects5 = ax.bar(ind - 0.27 + 2*width, np.array(agg_unique['kings_mov']), width, label='3-kings-mov',color='#058E41')
        rects7 = ax.bar(ind - 0.27 + 3*width, np.array(agg_unique['tc_mov']), width, label='TC-mov',color='#E5AB14')
        rects9 = ax.bar(ind - 0.27 + 4*width, np.array(agg_unique['cp_score']), width, label='CO-score',color='#9E6CC2')
        
        # Add some text for labels, title and custom x-axis tick labels, etc.
        #ax.set_ylabel('Average Number')
        ax.set_title(title[model],loc='right',size=14,fontname='Times')
        ax.set(ylim=[0,21])
        ax.set_xticks(ind)
        ax.set_xticklabels(ind_str)
        ax.tick_params(labelsize=13)
        ax.minorticks_on()
        ax.legend(loc='upper left')
        plt.savefig('unique-' + model + '2',dpi=500)
        plt.close()
