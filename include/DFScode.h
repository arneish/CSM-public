/* 
 * File:   DFScode.h
 * Author: Arneish Prateek
 * Created on 19 February, 2018, 10:19 PM
 */

#ifndef DFSCODE_H
#define DFSCODE_H
#include <vector>
#include <unordered_map>
#include <iostream>
#include <set>
#include <string>

using namespace std;

typedef vector<vector<pair<int, int>>> ADJ_LIST;

struct DFSCodeNode
{
    int a, b;
    int la, lab, lb;
    int vid_a, vid_b;

    DFSCodeNode(int _a = -1, int _b = -1, int _la = -1, int _lab = -1, int _lb = -1, int _vid_a = -1, int _vid_b = -1) : a(_a), b(_b), la(_la), lab(_lab), lb(_lb), vid_a(_vid_a), vid_b(_vid_b) {}
    ~DFSCodeNode() {}
    bool isForward() const
    {
        return a < b;
    }
    bool isBackward() const
    {
        return a > b;
    }
    bool operator<(const DFSCodeNode &o) const
    {
        if (this->isBackward() && o.isForward())
            return 1;
        else if (this->isBackward() && o.isBackward() && b < o.b)
            return 1;
        else if (this->isBackward() && o.isBackward() && b == o.b && lab < o.lab)
            return 1;
        else if (this->isForward() && o.isForward() && a > o.a)
            return 1;
        else if (this->isForward() && o.isForward() && a == o.a && la < o.la)
            return 1;
        else if (this->isForward() && o.isForward() && a == o.a && la == o.la && lab < o.lab)
            return 1;
        else if (this->isForward() && o.isForward() && a == o.a && la == o.la && lab == o.lab && lb < o.lb)
            return 1;
        return 0;
    }
    bool operator==(const DFSCodeNode &o) const
    {
        return a == o.a && b == o.b && la == o.la && lab == o.lab && lb == o.lb;
    }
    string show()
    {
        string s = to_string(a) + ',' + to_string(b) + ',' + to_string(la) + ',' + to_string(lab) + ',' + to_string(lb) + '.';
        return s;
    }
};

class DFSCode
{
  public:
    vector<DFSCodeNode> dfslist;
    unordered_map<int, int> visited;
    string dfscode_str;

  public:
    DFSCode *GenMin(int vnum, vector<vector<pair<int, int>>> &graph, int seq, unordered_map<int, int> track, vector<int> &vtx2lbl);
    DFSCode *GlobalMin(vector<vector<pair<int, int>>> &, vector<int> &v);
    vector<int> right_path();
    void DFSCodeSetString();
    void append(DFSCode *);
    DFSCode MinDFS(DFSCode[], int);
    int reorder(int, int);
    bool operator<(const DFSCode &o) const
    {
        int minsize = min(dfslist.size(), o.dfslist.size());
        for (int i = 0; i < minsize; i++)
        {
            if (dfslist[i] < o.dfslist[i])
                return 1;
            else if (dfslist[i] == o.dfslist[i])
                continue;
            else
                return 0;
        }
        return dfslist.size() > o.dfslist.size();
    }
    bool operator==(const DFSCode &o) const
    {
        if (dfslist.size() != o.dfslist.size())
            return 0;
        for (int i = 0; i < (int)dfslist.size(); i++)
            if (!(dfslist[i] == o.dfslist[i]))
                return 0;
        return 1;
    }
};

#endif /* DFSCODE_H */
