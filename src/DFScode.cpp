#include <vector>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstdio>
#include <ctime>
#include "DFScode.h"

using namespace std;

void DFSCode::DFSCodeSetString()
{
    this->dfscode_str = "";
    for (auto &i : dfslist)
    {
        this->dfscode_str += i.show();
    }
}

void DFSCode::append(DFSCode *D) //appending the argument's DFS code into the dfslist of the current class object
{
    if (D == nullptr)
    {
        return;
    }
    for (auto &i : D->dfslist)
    {
        this->dfslist.push_back(i);
    }
}

int DFSCode::reorder(int lastindex, int rootid)
{
    //input is the last index of forward edge of previous *i
    //returns the last index in current *i
    vector<DFSCodeNode> clist = this->dfslist;
    int count = 0;
    vector<DFSCodeNode> V;
    for (int j = 0; j < clist.size(); j++)
    {
        DFSCodeNode i = clist[j];
        if (count == 0)
        {
            DFSCodeNode *temp = new DFSCodeNode(rootid, i.b + lastindex - rootid, i.la, i.lab, i.lb, i.vid_a, i.vid_b);
            V.push_back(*temp);
            count = temp->b;
        }
        else
        {
            DFSCodeNode *temp = new DFSCodeNode(i.a + lastindex - rootid, i.b + lastindex - rootid, i.la, i.lab, i.lb, i.vid_a, i.vid_b);
            V.push_back(*temp);
            count = temp->b;
        }
    }
    this->dfslist = V;
    return count;
}

DFSCode DFSCode::MinDFS(DFSCode D[], int size)
{
    DFSCode ret = D[0];
    int len = D[0].dfslist.size();
    for (int i = 1; i < size; i++)
    {
        for (int j = 0; j < len; j++)
        {
            if (D[i].dfslist.at(j) < ret.dfslist.at(j))
            {
                ret = D[i];
                break;
            }
        }
    }
    return ret;
}

DFSCode* DFSCode::GlobalMin(vector<vector<pair<int, int>>> &graph, vector<int> &vertices_to_labels)
{
    /*returns minimum DFScode (string)*/
    vector<DFSCode *> candidate;
    int min_label = vertices_to_labels[0];
    for (int i = 0; i < graph.size(); i++)
    {
        if (vertices_to_labels[i] < min_label)
            min_label = vertices_to_labels[i];
    }
    vector<int> minlabelvtxID;
    for (int i = 0; i < graph.size(); i++)
    {
        if (vertices_to_labels[i] == min_label)
            minlabelvtxID.push_back(i);
    }
    for (auto &i : minlabelvtxID)
    {
        DFSCode *trial = new DFSCode();
        unordered_map<int, int> t;
        candidate.push_back(trial->GenMin(i, graph, 0, t, vertices_to_labels));
    }
    sort(candidate.begin(), candidate.end(), [](const DFSCode *pa, const DFSCode *pb) { return (*pa) < (*pb); });
    candidate[0]->DFSCodeSetString();
    for (int i=1 ; i<candidate.size(); i++)
    {
        delete candidate[i];
    }
    return candidate[0];
}

vector<int> DFSCode::right_path()
{
    vector<DFSCodeNode> &v = this->dfslist;
    int counter = -1;
    int search = 0, idtoadd = -1;
    vector<int> result;
    int j = 0;
    while (counter != v.size() - 1)
    {
        for (j = v.size() - 1; j >= counter; j--)
        {
            if (v[j].a == search)
            {
                counter = j;
                search = v[j].b;
                idtoadd = v[j].vid_a;
                break;
            }
        }
        result.push_back(idtoadd);
    }
    result.push_back(v[j].vid_b);
    reverse(result.begin(), result.end()); //Rightmost path (vertex ID's in priority order)
    return result;
}

DFSCode *DFSCode::GenMin(int vnum, vector<vector<pair<int, int>>> &graph, int seq, unordered_map<int, int> track, vector<int> &vtx2lbl)
{
    track[vnum] = 1;
    vector<DFSCodeNode *> vecf;
    vector<DFSCode *> candidate;
    vector<int> reset_track;
    int fsize = 0;
    int basecase = -1;
    for (auto &i : graph[vnum])
    {
        if (track.count(i.first) == 0)
        {
            basecase = 0;
            reset_track.push_back(i.first);
            DFSCodeNode *temp = new DFSCodeNode(seq, seq + 1, vtx2lbl[vnum], i.second, vtx2lbl[i.first], vnum, i.first);
            DFSCode *subtree = new DFSCode();
            subtree->dfslist.push_back(*temp);
            subtree->append(GenMin(i.first, graph, seq + 1, track, vtx2lbl));
            candidate.push_back(subtree);
        }
    }
    if (basecase == -1)
    {
        return nullptr;
    }
    for (auto &i : reset_track)
    {
        track.erase(i);
    }
    reset_track.clear();
    this->dfslist.clear();
    sort(candidate.begin(), candidate.end(), [](const DFSCode *pa, const DFSCode *pb) { return (*pa) < (*pb); });
    int lastindex = seq;
    for (auto &i : candidate)
    {
        //Updating DFS Lists of candidates in sorted order using reorder():
        lastindex = i->reorder(lastindex, seq);
        this->append(i);
    }
    return this;
}
