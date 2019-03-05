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
#include <string>
#include <queue>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include <assert.h>
#include <unistd.h>

#include "DFScode-feb-22.h"

using namespace std;

#define LL long long int
#define PAUSE_TIME 3

struct edge;
struct replica_s;
struct subgraph_pattern;
struct correlated_patterns;
class ComparePattern;
class CompareCorrelation;

extern ofstream frequent_txt;
extern ofstream pfile;

extern unordered_map<int, int> patternSizeToCount;
extern unordered_map<int, int> patternSizeToCount2;
extern unsigned int pause_;              /*To invoke sleep() for visualising values during debugging*/
extern vector<int> pattern_sizes;        /*To store pattern sizes of explored patterns*/
extern unordered_set<string> dictionary; /*To store DFS labelling of explored graphs*/
extern vector<set<int> *> corvm;         /*To store corV set for every vertex*/
extern vector<set<int>> corP;            /*To store corP set for every vertex*/
extern vector<int> vertices_to_labels;   /*To store mappings from vertex to label for input graph*/
extern int global_num_patterns;          /*To store number of patterns explored and assign pattern ID*/
extern vector<bool> frequent_vertices;   /*To store whether a vertex is frequent or not*/
extern int k;                            /*'K' in Top-K*/
extern int min_sup;                      /*Support threshold*/
extern int control;                      /*To control complete or approximate execution*/
extern int hop;                          /*The 'hop' parameter*/
extern vector<int> freq_vert_patternind;
extern int global_dict_hits;
extern vector<vector<int>> label_to_vertices;
/*To store vertex IDs for every label in the input graph*/
extern vector<vector<pair<int, int>>> inp_graph;
/*To store the input graph, pair<a,b> represents edge to vertex-ID 'a' with edge-label 'b'*/

extern priority_queue<subgraph_pattern, vector<subgraph_pattern>, ComparePattern> leaf;
extern vector<subgraph_pattern> operated;
extern priority_queue<correlated_patterns, vector<correlated_patterns>, CompareCorrelation> top_k;

extern void freq_vertices();
extern inline bool NotFrequent(pair<int, int> &v_pair)
{
    return (frequent_vertices[v_pair.first] == false);
}
extern void modify_adjlist();
extern void compute_corv(vector<vector<pair<int, int>>> &, vector<set<int> *> &, int &);
extern void corvtraverse(vector<vector<pair<int, int>>> &, vector<set<int> *> *, vector<set<int> *> *);
extern void search();
extern inline bool ceasing_condition();
extern void operate(subgraph_pattern &);
extern void constructcollectstat(subgraph_pattern &, int, unordered_map<int, set<int>> &);
extern void allinstances_collectstat(subgraph_pattern &, vector<struct edge *> *, int &, int, unordered_set<int> &);
extern int DFS_instance_all_collectstat(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, unordered_set<int> &);
extern void extend(subgraph_pattern &);
extern void possible_extensions_(subgraph_pattern &, vector<LL> &);
extern bool instance_checkerDFS(subgraph_pattern &, vector<struct edge *> *, int &, const int &, int &);
extern bool DFS_instance(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &instance, set<int> &, int &);
/*instance-map: (pID, replica-ID) constituting that instance*/
extern pair<int, int> support_index(replica_s *);
/*Set-A CSM-approximate replica-extension:*/
extern replica_s *get_replica_brute(subgraph_pattern &, int, int, int);
extern void instance_allDFS(subgraph_pattern &, vector<struct edge *> *, int &, const int &, const int &, vector<unordered_set<int>> &, bool &);
extern int DFS_instance_all(subgraph_pattern &, int, vector<struct edge *> *, unordered_map<int, int> &, set<int> &, const int &, vector<unordered_set<int>> &, bool &);
/*Set-B CSM-complete replica-extension:*/
extern replica_s *get_replica_brute_orig(subgraph_pattern &q, int extending_index, int extended_label, int edge_label);
extern void instance_allDFS_orig(subgraph_pattern &q, vector<struct edge *> *dfslist, int &extending_index, int extending_index_vertexID, int &conflict_ID, vector<unordered_set<int>> &new_mappings, bool &found);
extern void DFS_instance_all_orig(subgraph_pattern &q, int list_index, vector<struct edge *> *dfslist, unordered_map<int, int> &instance, set<int> &instance_set, int &conflict_ID, vector<unordered_set<int>> &new_mappings, bool &found);
extern void print_pattern (vector<vector<pair<int, int>>> &, vector<int> &);
extern void write_pattern(vector<vector<pair<int, int>>> &, vector<int> &, ofstream &, int, string);

struct edge /**/
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
    unordered_map<string, unordered_set<int>> first_visit;             //key=pattern-vertex ID, value={graph vertex-IDs visited once during get_replica_brute exploration}
    unordered_map<string, unordered_set<int>> potential_child_mapping; //key=graph-vertex ID, value={{child-pattern-vertex ID,{potentially valid graph-vertex IDs}}}
    unordered_map<string, unordered_set<int>> confirmed_child_mapping;
    //unordered_map<string, unordered_set<int>> remove_potential_child_mapping; // for_updation of previous map
    vector<vector<int>> parent_to_child;
    vector<int> leaf_vertexID;
    vector<unordered_set<int>> enumerated_mappings;

    unordered_map<int, set<int>> centerToVertexSet; //key=center(graph-vertex)ID, value={vertex-IDs constituting instance group} [introduced Nov 14, 2018]
    unordered_map<int, int> centerToClusterCenter;  //key=center(graph-vertex)ID, value=cluster center(graph-vertex)ID: Initialised as self in clustercollectstat
    unordered_map<int, set<int>> vertexToGroups;    //key=graph-vertex ID, value={Instance group center(graph-vertex) IDs}

    vector<unordered_set<int>> mappings;
    unordered_map<int, vector<pair<int, int>>> adj_list; //key=graph-vertex ID, value={pair(graph-vertex ID, edge-label)}
    unordered_map<int, set<int>> vmaplist;               //key=vertex ID, value={pattern vertex IDs}: Modifications on June 10, 2018 by A.P.

    replica_s(unordered_map<int, vector<pair<int, int>>> &adj_list_, int pattern_size)
    {
        this->adj_list = adj_list_;
    }


    replica_s(vector<unordered_set<int>> &mappings, unordered_map<int, vector<pair<int, int>>> &adj_list_, unordered_map<int, set<int>> &vmaplist)
    {
        this->mappings = mappings;
        this->adj_list = adj_list_;
        this->vmaplist = vmaplist;
    }
};

struct subgraph_pattern /*The structure for a subgraph pattern*/
{
    int id;
    vector<vector<pair<int, int>>> graph;
    vector<int> vertices_to_labels;
    replica_s *replica;
    int support;
    vector<int> ordering;
    int min_index;
    vector<int> rightmost_path;
    string strdfscode;
    set<int> subgraphs;
    //Modification: August 15: BFSorder storage for different roots (key), value (vector of BFS_order: [element, parent])-on-demand:
    unordered_map<int, vector<pair<int, int>>> graph_bfsorder;
    unordered_map<int, vector<struct edge *>> graph_dfsorder; //key: root, value: vector of DFS_order <childID, parentID> with parentID of root = -1

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

    //modification Aug 15: BFS_order construction for on-demand calling & return:
    vector<pair<int, int>> bfsordercons(int id)
    {
        if (this->graph_bfsorder.count(id))
            return this->graph_bfsorder[id];
        vector<pair<int, int>> order = {make_pair(id, -1)}; //initialised with root, parent(=-1)
        unordered_map<int, int> explored_set;
        int index = 0;
        while (index < order.size())
        {
            for (auto &elem : graph[order[index].first])
            {
                if (!explored_set.count(elem.first))
                {
                    order.push_back(make_pair(elem.first, order[index].first));
                    explored_set[elem.first] = 1;
                }
            }
            index++;
        }
        this->graph_bfsorder[id] = order;
        return this->graph_bfsorder[id];
    }

    //DFS_order construction for on-demand calling and return;
    vector<struct edge *> *dfsordercons(int id)
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
