/* 
 * File:   csm.h
 * Author: Arneish Prateek
 *
 * Created on 21 February, 2019, 8:53 PM
 */
#ifndef CSM_H
#define CSM_H
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <string>
#include <queue>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include <stdio.h>
#include <cstdio>
#include <assert.h>
#include <unistd.h>
#include "DFScode.h"

using namespace std;

#define LL long long int
typedef pair<int, int> edge_tuple;

struct edge;
struct replica_s;
struct subgraph_pattern;
struct correlated_patterns;
class ComparePattern;
class CompareCorrelation;

extern clock_t project_start_time;
extern int search_counter;
extern ofstream output_txt;                        /*To store top-K results and operated patterns*/
extern unordered_map<int, int> patternSizeToCount; /*To store count of different pattern sizes explored*/
extern unsigned int pause_;                        /*To invoke sleep() for visual debugging*/
extern unordered_set<string> dictionary;           /*To store DFS labelling of explored graphs*/
extern vector<set<int> *> corvm;                   /*To store corV set for every vertex*/
extern vector<set<int>> corP;                      /*To store corP set for every vertex*/
extern vector<int> vertices_to_labels;             /*To store mappings from vertex to label for input graph*/
extern int global_num_patterns;                    /*To store number of patterns explored and assign pattern ID*/
extern vector<bool> frequent_vertices;             /*To store whether a vertex is frequent or not*/
extern int k;                                      /*'K' in Top-K*/
extern int min_sup;                                /*Support threshold*/
extern int control;                                /*To control CSM-complete or CSM-approximate execution*/
extern int hop;                                    /*The 'hop' parameter*/
extern vector<int> freq_vert_patternind;           /*To store patterns ID's of every frequent vertex*/
extern int global_dict_hits;                       /*To store number of patterns that matched a dictionary pattern*/
extern vector<vector<int>> label_to_vertices;      /*To store vertex IDs for every label in the input graph*/
extern vector<vector<pair<int, int>>> inp_graph;
/*To store the input graph, pair<a,b> represents edge to vertex-ID 'a' with edge-label 'b'*/
extern priority_queue<subgraph_pattern, vector<subgraph_pattern>, ComparePattern> leaf;
/*To store frequent patterns in a priority queue for subsequent "operation", Comparator: "ComparePattern"*/
extern vector<subgraph_pattern> operated;
/*To store frequent subgraph patterns that have been "operated"*/
extern priority_queue<correlated_patterns, vector<correlated_patterns>, CompareCorrelation> top_k;
/*To store a priority queue of top-K correlated patterns, Comparator: "CompareCorrelation"*/
extern vector<vector<int>> gs_instances;

/*Main CSM functions:*/
extern void freq_vertices();
extern inline bool NotFrequent(pair<int, int> &v_pair)
{
    return (frequent_vertices[v_pair.first] == false);
}
extern void modify_adjlist();
extern void compute_corv(vector<vector<pair<int, int>>> &, vector<set<int> *> &, int &);
extern void corvtraverse(int, int, vector<vector<pair<int, int>>> &, vector<set<int>*> *, vector<set<int>*> *);
extern void search();
extern inline bool ceasing_condition();
extern void operate(subgraph_pattern &);

/*Collect-stat functions:*/
/*Set-A CSM-approximate*/
extern void constructcollectstat(subgraph_pattern &, int, unordered_map<int, set<int>> &);
extern void allinstances_collectstat(subgraph_pattern &, vector<struct edge *> *, int &, int, unordered_set<int> &);
extern bool DFS_instance_all_collectstat(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, unordered_set<int> &);
/*Set-B CSM-complete*/
extern void constructcollectstat_orig(subgraph_pattern &, int, unordered_map<int, set<int>> &);
extern void allinstances_collectstat_orig(subgraph_pattern &, vector<struct edge *> *, int &, int, unordered_set<int> &);
extern void DFS_instance_all_collectstat_orig(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, unordered_set<int> &);

/*Pattern extension functions:*/
extern void extend(subgraph_pattern &);
extern void possible_extensions_(subgraph_pattern &, unordered_set<LL> &);//vector<LL> &);
extern bool instance_checkerDFS(subgraph_pattern &, vector<struct edge *> *, int &, const int &, int &);
extern bool DFS_instance(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, int &);
/*instance-map: (pID, replica-ID) constituting that instance*/
extern pair<int, int> support_index(replica_s *);

/*Set-A CSM-approximate replica-extension:*/
extern replica_s *get_replica_brute(subgraph_pattern &, int, int, int);
extern bool instance_allDFS(subgraph_pattern &, vector<struct edge *> *, int &, const int &, const int &, vector<unordered_set<int>> &);
extern bool DFS_instance_all(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, const int &, vector<unordered_set<int>> &);
/*Set-B CSM-complete replica-extension:*/
extern replica_s *get_replica_brute_orig(subgraph_pattern &, int, int, int);
extern void instance_allDFS_orig(subgraph_pattern &, vector<struct edge *> *, int &, int, int &, vector<unordered_set<int>> &, bool &);
extern void DFS_instance_all_orig(subgraph_pattern &, int list_index, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, int &, vector<unordered_set<int>> &, bool &);
extern void print_pattern(vector<vector<pair<int, int>>> &, vector<int> &);
extern void write_pattern(vector<vector<pair<int, int>>> &, vector<int> &, ofstream &, int, string);

struct edge /*To be used in the DFSList vector to store pattern edges*/
{
    int parent;
    int child;
    int edge_label;

    edge(int &parent, int &child, int &edge_label)
    {
        this->parent = parent;
        this->child = child;
        this->edge_label = edge_label;
    }
};

struct replica_s /*The replica structure*/
{
    unordered_map<string, unordered_set<int>> first_visit;             /*key=pattern vertex ID, value={graph vertex IDs visited once during get_replica_brute exploration}*/
    unordered_map<string, unordered_set<int>> potential_child_mapping; /*key={child vertex ID}, value={potential graph vertex IDs, mapping to the child vertex}*/
    vector<unordered_set<int>> enumerated_mappings;
    /*To store graph vertices mapping to pattern vertex at the corresponding ID such that all mappings of its child are enumerated*/
    unordered_map<string, unordered_set<int>> confirmed_child_mapping; /*key={child vertex ID}, value={confirmed graph vertex IDs, mapping to the child vertex}*/
    vector<vector<int>> parent_to_child;                               /*To store children of parent at corresponding index in the DFS tree*/
    vector<int> leaf_vertexID;                                         /*To store patterns at the leaves of the DFS tree*/
    vector<unordered_set<int>> mappings;                               /*To store graph vertices that map to pattern vertices at corresponding index*/
    unordered_map<int, set<edge_tuple>> adj_list;                      /*key=graph vertex ID, value={pair(graph vertex ID, edge-label)}*/
    unordered_map<int, set<int>> vmaplist;                             /*key=graph vertex ID, value={pattern vertex IDs}*/

    replica_s(unordered_map<int, set<edge_tuple>> &adj_list_)
    {
        this->adj_list = adj_list_;
    }

    replica_s(vector<unordered_set<int>> &mappings, unordered_map<int, set<edge_tuple>> &adj_list_, unordered_map<int, set<int>> &vmaplist)
    {
        this->mappings = mappings;
        this->adj_list = adj_list_;
        this->vmaplist = vmaplist;
    }
};

struct subgraph_pattern /*The structure for a subgraph pattern*/
{
    int id;                                                   /*Pattern IDa*/
    vector<vector<pair<int, int>>> graph;                     /*To store pattern graph*/
    vector<int> vertices_to_labels;                           /*To store labels for the vertex at the corresponding index*/
    replica_s *replica;                                       /*Mappings of this->graph in the input graph*/
    int support;                                              /*MNI-support for this this->graph*/
    int min_index;                                            /*Vertex with minimum number of distinct mappings*/
    string strdfscode;                                        /*DFS code for this->graph*/
    vector<int> rightmost_path;                               /*Rightmost-path vertices for this->graph*/
    set<int> subgraphs;                                       /*Set of IDs of subgraphs of this->graph*/
    unordered_map<int, vector<struct edge *>> graph_dfsorder; /*key=root, value=vector of DFS_order <childID, parentID> with parentID of root = -1*/

    subgraph_pattern(int &id, vector<vector<pair<int, int>>> &graph_, vector<int> &vertices_to_labels, replica_s *rep, int &sup, int &min_index, string &dfscode_, vector<int> &rm_path)
    {
        this->id = id;
        this->graph = graph_;
        this->vertices_to_labels = vertices_to_labels;
        this->replica = rep;
        this->support = sup;
        this->min_index = min_index;
        this->strdfscode = dfscode_;
        this->rightmost_path = rm_path;
    }

    vector<struct edge *> *dfsordercons(int id) /*DFS_order construction for on-demand calling and return*/
    {
        if (this->graph_dfsorder.count(id))
            return &this->graph_dfsorder[id];
        vector<struct edge *> dfslist;
        unordered_set<int> explored_set;
        DFS(id, -1, -1, dfslist, explored_set);
        this->graph_dfsorder[id] = dfslist;
        return &this->graph_dfsorder[id];
    }

    void DFS(int id, int parentID, int edge_label, vector<struct edge *> &dfslist, unordered_set<int> &explored_set)
    {
        struct edge *e = new edge(parentID, id, edge_label);
        dfslist.push_back(e);
        explored_set.insert(id);
        for (auto &children : this->graph[id])
        {
            if (!explored_set.count(children.first))
                DFS(children.first, id, children.second, dfslist, explored_set);
        }
    }
};

struct correlated_patterns
{
    int operated_id_1;
    int operated_id_2;
    int corr_value;
    correlated_patterns(int &t1, int &t2, int &t3)
    {
        this->operated_id_1 = t1;
        this->operated_id_2 = t2;
        this->corr_value = t3;
    }
};

class ComparePattern
{
  public:
    bool operator()(subgraph_pattern &t1, subgraph_pattern &t2) // Returns true if t1 is earlier than t2
    {
        if (t1.support < t2.support)
            return true;
        else if (t1.support == t2.support)
        {
            if (t1.replica->mappings.size() > t2.replica->mappings.size())
                return true;
            return false;
        }
        return false;
    }
};

class CompareCorrelation
{
  public:
    bool operator()(correlated_patterns &t1, correlated_patterns &t2)
    {
        if (t1.corr_value > t2.corr_value)
            return true;
        return false;
    }
};

#endif /* CSM_H */
