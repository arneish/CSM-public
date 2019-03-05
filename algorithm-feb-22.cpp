#include "csm.h"

using namespace std;

void freq_vertices()
{
	for (int i = 0; i < inp_graph.size(); i++)
	{
		frequent_vertices.push_back(false); //Initially, all vertices are NOT frequent
	}
	cerr << "Label_to_vertices size:" << label_to_vertices.size() << "\n";
	vector<int> label_to_freq(static_cast<int>(label_to_vertices.size()), 0);
	int max_mapping = 0;
	for (int i = 0; i < label_to_vertices.size(); i++) //labels being identified as [0, label_to_vertices.size()-1]
	{
		max_mapping = max(max_mapping, static_cast<int>(label_to_vertices[i].size()));
		label_to_freq[i] = label_to_vertices[i].size();
		if (label_to_freq[i] >= min_sup)
		{
			cout << "Frequent vertex-label: " << i << " with frequency: " << label_to_vertices[i].size() << "\n";
			int id = global_num_patterns++;
			freq_vert_patternind.push_back(id);

			/*Constructing and storing the replica for the frequent pattern:*/
			unordered_set<int> maps;
			unordered_map<int, vector<pair<int, int>>> replica_adj_list;
			unordered_map<int, set<int>> vmaplist;
			for (int j = 0; j < label_to_vertices[i].size(); j++) //iterating over all vertices of the frequent label[i]
			{
				frequent_vertices[label_to_vertices[i][j]] = true;
				set<int> v_map = {0}; //pattern-vtxID is trivially zero
				vmaplist[label_to_vertices[i][j]] = v_map;
				vector<pair<int, int>> empty_adj_list_;
				replica_adj_list[label_to_vertices[i][j]] = empty_adj_list_;
				maps.insert(label_to_vertices[i][j]); // "maps": stores all vtxID's of freq label[i]
			}
			vector<unordered_set<int>> mappings = {maps};
			replica_s *replica_ = new replica_s(mappings, replica_adj_list, vmaplist); //"replica_": replica for the frequent label

			/*Constructing the frequent pattern graph and the subgraph_pattern data structure:*/
			vector<pair<int, int>> empty_adj_list;
			vector<vector<pair<int, int>>> graph_ = {empty_adj_list};
			vector<int> vertices_to_labels = {i};
			int support = label_to_vertices[i].size();
			int min_index = 0;
			string dfscode_ = to_string(i); //DFS-code of vertices is the vertex label
			vector<int> rightmost_path = {0};
			subgraph_pattern pattern_(id, graph_, vertices_to_labels, replica_, support, min_index, dfscode_, rightmost_path);
			leaf.push(pattern_);
		}
	}
	sort(label_to_freq.begin(), label_to_freq.end(), greater<int>());
	int counter = 0;
	cerr << "\nLabel Index: "
		 << "Frequency\n";
	while (counter < label_to_freq.size())
	{
		if (label_to_freq[counter] >= min_sup)
			cerr << counter << ": " << label_to_freq[counter] << " ";
		counter++;
	}
	cerr << "\nInitial leaf size: " << leaf.size() << "\n";
	cerr << "Max mappings size: " << max_mapping << "\n";
	cerr << "Returning from freq_vertices()\n";
}

void modify_adjlist()
{
	//NOTE: only deletion from Adjacency Lists; no vertex deletion to maintain original vertex-ID's
	int edges_deleted = 0, pre_size, post_size;
	for (auto &i : inp_graph)
	{
		pre_size = i.size();
		i.erase(remove_if(i.begin(), i.end(), NotFrequent), i.end());
		post_size = i.size();
		edges_deleted += (pre_size - post_size);
	}
	cerr << "edges_deleted = " << edges_deleted << "\n";
}

void compute_corv(vector<vector<pair<int, int>>> &graph, vector<set<int> *> &corv, int &hop)
{
	/*//hop = 1 code:
    if (hop == 1) 
    {
        for (int i = 0; i < graph.size(); i++)
        {
            set<int> *temp = new set<int>;
            for (auto &j : graph[i])
            {
                temp->insert(j.first);
            }
            corv.push_back(temp);
        }
        return;
    } */
	/*GENERAL CASE: hop >= 1 code:*/
	vector<vector<set<int> *> *> list;
	/*The i-th vector element is a hashmap of vertex-id to id's of its L-<i-th> neighbours
      L-0 ("Level"-0) for a vertex is the vertex itself, L-1 is the set of immediate neighbours + itself and so on. */
	list.push_back(new vector<set<int> *>);
	for (int i = 0; i < graph.size(); i++)
	{
		set<int> *L1 = new set<int>{i};
		for (auto &j : graph[i])
		{
			L1->insert(j.first);
		}
		list[0]->push_back(L1); //L-1 neighbours stored in list
	}
	for (int i = 1; i < hop; i++)
	{
		cerr << "hop#: " << i + 1 << "\n";
		list.push_back(new vector<set<int> *>);
		corvtraverse(graph, list[i - 1], list[i]);
		cerr << "clearing list: L-" << i << "\n";
		for (auto &corv_set : *list[i - 1])
		{
			delete corv_set;
		}
		list[i - 1]->clear();
		delete list[i - 1];
	}
	int total_corv_size = 0;
	for (int i = 0; i < graph.size(); i++) //inserting into corv:
	{
		set<int> *corv_members = list[hop - 1]->at(i);
		corv.push_back(corv_members);
		total_corv_size += corv_members->size();
	}
	cerr << "Total corV size: " << total_corv_size << "\n";
	// cerr << "UNIX Pause_\n";
	// sleep(pause_);
	// sleep(pause_);
}

void corvtraverse(vector<vector<pair<int, int>>> &graph, vector<set<int> *> *corvmap, vector<set<int> *> *corvmap2)
{
	for (int i = 0; i < graph.size(); i++)
	{
		cerr << "\r" << i << "/" << graph.size() << " ";
		corvmap2->push_back(new set<int>{i});
		for (auto &nbr : graph[i])
		{
			corvmap2->at(i)->insert(corvmap->at(nbr.first)->begin(), corvmap->at(nbr.first)->end());
		}
	}
	cerr << "\n";
}

void search()
{
	int count = 0;
	while (!ceasing_condition())
	{
		cerr << "Ceasing condition not met. Search loop counter: " << ++count
			 << "\nLeaf Size: " << leaf.size() << ", Top-k size: " << top_k.size() << "/" << k << "\n";
		subgraph_pattern q = leaf.top();
		leaf.pop();
		patternSizeToCount[q.graph.size()] += 1;
		operate(q);
		cerr << "Post operate() Top-k size: " << top_k.size() << "\n";
		extend(q);
		cerr << "Post extend() Leaf size:" << leaf.size() << "\n";
	}
}

inline bool ceasing_condition()
{
	if (top_k.size() == k && leaf.top().support <= top_k.top().corr_value)
		return true;
	else if (leaf.empty())
		return true;
	return false;
}

void operate(subgraph_pattern &q)
{
	int operated_size = operated.size();
	cerr << "operate() call:\n";
	unordered_map<int, set<int>> collection;
	//"collection": key: ID of mapping of center-ID, value: set of proximal patternIDs for that mapping
	constructcollectstat(q, q.min_index, collection);
	for (int i = 0; i < operated_size; i++)
	{
		if (q.subgraphs.find(operated[i].id) == q.subgraphs.end() && q.replica->mappings.size() > 1)
		{
			int pattern_ind = operated[i].id;
			unordered_set<int> &root_mappings = q.replica->mappings[q.min_index];
			int count = 0;
			for (auto it : root_mappings)
			{
				if (collection[it].find(pattern_ind) != collection[it].end())
					count++;
			}
			if (top_k.size() < k && count > 0)
			{
				struct correlated_patterns *patterns_pair = new struct correlated_patterns(i, operated_size, count);
				top_k.push(*patterns_pair);
			}
			else if (top_k.top().corr_value < count)
			{
				top_k.pop();
				struct correlated_patterns *patterns_pair = new struct correlated_patterns(i, operated_size, count);
				top_k.push(*patterns_pair);
			}
		}
	}
	operated.push_back(q);
}

void constructcollectstat(subgraph_pattern &q, int centerID, unordered_map<int, set<int>> &collectSetMap)
{
	vector<vector<pair<int, int>>> &Q_graph = q.graph;
	vector<unordered_set<int>> &Q_mappings = q.replica->mappings;
	/*construct DFS-order from centerID*/
	vector<struct edge *> *dfslist = q.dfsordercons(centerID);
	q.replica->parent_to_child.clear();
	q.replica->parent_to_child.resize(q.graph.size(), vector<int>());
	for (int i = 1; i < dfslist->size(); i++)
	{
		q.replica->parent_to_child[dfslist->at(i)->parent].push_back(dfslist->at(i)->child);
	}
	vector<int> &leaves = q.replica->leaf_vertexID;
	leaves.clear();
	for (int i = 0; i < q.replica->parent_to_child.size(); i++)
	{
		if (!q.replica->parent_to_child[i].size())
		{
			leaves.push_back(i);
		}
	}
	q.replica->enumerated_mappings.clear();
	q.replica->enumerated_mappings.resize(q.graph.size(), unordered_set<int>());
	int map_counter = 0;
	unordered_set<int> group_members; /*to store all vertex-IDs constituting an instance group corresponding to an instance center*/
	for (const auto &mapped_centerID : Q_mappings[centerID])
	{
		group_members.clear();
		if (q.graph.size() >= 3)
		{
			cerr << "\r" << map_counter++ << "/" << Q_mappings[centerID].size();
		}
		allinstances_collectstat(q, dfslist, centerID, mapped_centerID, group_members);
		/*
		1. Aggregate proximal pattern IDs for every instance group
		2. Assign self-patternID to all vertices within the hop of every vertex constituting the instance group 
		*/
		set<int> &proximal_patterns = collectSetMap[mapped_centerID];
		for (const auto &vertex : group_members)
		{
			proximal_patterns.insert(corP[vertex].begin(), corP[vertex].end());
			for (auto &corv_vertex : *corvm[vertex])
			{ //put self-patternID in all of corV for this vertex
				corP[corv_vertex].insert(q.id);
			}
		}
	}
	/*
	(Optional)Reset q.replica: first visit map, potential child map, confirmed child map 
	*/
	q.replica->first_visit.erase(q.replica->first_visit.begin(), q.replica->first_visit.end());
	q.replica->potential_child_mapping.erase(q.replica->potential_child_mapping.begin(), q.replica->potential_child_mapping.end());
	q.replica->confirmed_child_mapping.erase(q.replica->confirmed_child_mapping.begin(), q.replica->confirmed_child_mapping.end());
}

void allinstances_collectstat(subgraph_pattern &q, vector<struct edge *> *dfslist, int &center_ID, int center_vertexID, unordered_set<int> &group_members)
{
	unordered_map<int, int> instance; // <patternvtxID, replicavtxID (mapping)>
	instance[center_ID] = center_vertexID;
	set<int> instance_set = {center_vertexID};
	int list_index = 1; //start DFS-wise instance enumeration from dfslist[1] (dfslist[0] is <center-vertex, -1>)
	DFS_instance_all_collectstat(q, list_index, dfslist, instance, instance_set, group_members);
}

int DFS_instance_all_collectstat(subgraph_pattern &q, int list_index, vector<struct edge *> *dfslist, unordered_map<int, int> &instance, set<int> &instance_set, unordered_set<int> &group_members)
{
	/*returns true iff at least one instance found*/

	bool local_found = false; //"true" iff including current child_mapping results in at least one valid instance
	/* base case: */
	if (list_index == dfslist->size())
	{
		//dfslist has been traversed till the end
		for (auto &instance_mapping : instance)
		{
			group_members.insert(instance_mapping.second);
		}
		for (const auto &leaf_vertex : q.replica->leaf_vertexID)
		{
			q.replica->enumerated_mappings[leaf_vertex].insert(instance[leaf_vertex]);
		}
		return local_found = true;
	}
	/* recursion: */
	/*If first-visit then store set of all possible child-mapping vertices, else iterate over existing set of (remaining) potential child mappings*/
	int &parent_patternID = dfslist->at(list_index)->parent; //vertex-ID of parent vertex in the pattern graph
	int &parent_vertexID = instance[parent_patternID];
	int &child_patternID = dfslist->at(list_index)->child; //vertex-ID of child vertex (wrt DFS list order) in the pattern graph
	int &edge_label = dfslist->at(list_index)->edge_label;
	unordered_set<int> &potential_set = q.replica->potential_child_mapping[to_string(parent_vertexID) + ' ' + to_string(child_patternID)];
	unordered_set<int> &confirmed_set = q.replica->confirmed_child_mapping[to_string(parent_vertexID) + ' ' + to_string(child_patternID)];

	if (!q.replica->first_visit[to_string(parent_patternID) + ' ' + to_string(child_patternID)].count(parent_vertexID))
	{
		//this is the first-visit at pattern_vertexID for this edge-type
		for (const auto &edge : q.replica->adj_list[parent_vertexID])
		{
			if (q.replica->vmaplist[edge.first].count(child_patternID) && edge.second == edge_label)
			{
				potential_set.insert(edge.first);
				/*store the set of potentially-valid child mappings*/
			}
		}
		q.replica->first_visit[to_string(parent_patternID) + ' ' + to_string(child_patternID)].insert(parent_vertexID);
	}

	/*iterate over (remaining) potential child mappings and attempt to enumerate instances*/
	unordered_set<int> add_to_confirmed;
	for (auto &child_mapping : potential_set)
	{
		if (instance_set.count(child_mapping))
			continue;
		instance[child_patternID] = child_mapping;
		instance_set.insert(child_mapping);
		bool instanceFound = DFS_instance_all_collectstat(q, list_index + 1, dfslist, instance, instance_set, group_members);
		if (instanceFound)
			local_found = true;
		/*Check if child_mapping is "completely-enumerated"*/
		if (instanceFound && q.replica->enumerated_mappings[child_patternID].count(child_mapping))
		{
			add_to_confirmed.insert(child_mapping);
			//to be deleted from potential_child_mapping and move to confirmed_child_mapping
		}
		instance.erase(child_patternID);
		instance_set.erase(child_mapping);
	}

	if (!local_found)
	{
		/*in this case (i.e. not a single valid instance found using potential set), 
		iterate over confirmed child mappings and try to enumerate (one) instance*/
		for (auto &child_mapping : confirmed_set)
		{
			if (instance_set.count(child_mapping))
				continue;
			instance[child_patternID] = child_mapping;
			instance_set.insert(child_mapping);
			bool instanceFound = DFS_instance_all_collectstat(q, list_index + 1, dfslist, instance, instance_set, group_members);
			instance.erase(child_patternID);
			instance_set.erase(child_mapping);
			if (instanceFound)
			{
				local_found = true;
				break;
			}
		}
	}

	/*iterate over add_to_confirmed set and update potential_child_mapping (deletion) and confirmed_child_mapping (insertion) */
	confirmed_set.insert(add_to_confirmed.begin(), add_to_confirmed.end());
	for (auto &child_mapping : add_to_confirmed)
	{
		potential_set.erase(child_mapping);
	}
	add_to_confirmed.clear();

	/*if ALL left-over children sets are empty, parent-vertex is also completely enumerated FOR THE CORRESPONDING MAPPING */
	bool ToBeMarkedEnumerated = true;
	for (auto &child_pattern : q.replica->parent_to_child[parent_patternID])
	{
		string key = to_string(parent_vertexID) + ' ' + to_string(child_pattern);
		if (!q.replica->potential_child_mapping.count(key) || q.replica->potential_child_mapping[key].size())
		{
			ToBeMarkedEnumerated = false;
			break;
		}
	}
	if (ToBeMarkedEnumerated)
	{
		q.replica->enumerated_mappings[parent_patternID].insert(parent_vertexID);
	}
	return local_found;
}

void extend(subgraph_pattern &q)
{
	//try to call minDFS once only here
	cerr << "extend() call:\n";
	vector<vector<pair<int, int>>> &Q_graph = q.graph;
	vector<int> &Q_vertices_to_labels = q.vertices_to_labels;
	vector<LL> extensions;
	possible_extensions_(q, extensions);
	int extensions_size = extensions.size();
	cerr << "Candidate extensions size: " << extensions_size << "\n";
	pair<int, int> sup_ind_org = support_index(q.replica);
	int &support_org = sup_ind_org.first;
	int &min_index_org = sup_ind_org.second;
	for (int i = 0; i < extensions.size(); i++)
	{
		cerr << "Candidate extn #" << i + 1 << ": " << extensions[i] << "\n";
	}

	/*check size of q and same label or not*/
	// bool allLabelSame = true;
	// int label = q.vertices_to_labels[0];
	// for (int i = 1; i < q.graph.size(); i++)
	// {
	// 	if (label != q.vertices_to_labels[i])
	// 	{
	// 		allLabelSame = false;
	// 		break;
	// 	}
	// }

	for (int i = 0; i < extensions_size; i++)
	{
		cerr << "extension: " << i + 1 << "/" << extensions_size << ": ";
		LL extension_ = extensions[i];
		int vert_id = extension_ % 1000;
		extension_ /= 1000;
		int edge_label = extension_ % 10000;
		extension_ /= 10000;
		int extended_label = extension_;
		cerr << " vert_id: " << vert_id << " edge_l: " << edge_label << " extended_l: " << extended_label << "\n";
		// if (q.graph.size() >= 5 && allLabelSame)
		// {
		// 	if (extended_label == label)
		// 	{
		// 		cerr << "Extended_label same as all vertex labels - SKIPPING EXTENSION \n";
		// 		continue;
		// 	}
		// }
		// else if (q.graph.size() == 3)
		// {
		// 	continue;
		// }

		// if (q.graph.size() == 5)
		// {
		// 	continue;
		// }
		// if (q.graph.size() == 4 && allLabelSame)
		// 	continue;

		/* Construct the candidate pattern graph (tree): */
		vector<vector<pair<int, int>>> candidate_graph = q.graph;
		int N = candidate_graph.size();
		candidate_graph[vert_id].push_back(make_pair(N, edge_label));
		vector<pair<int, int>> new_vertex_list = {make_pair(vert_id, edge_label)};
		candidate_graph.push_back(new_vertex_list);
		vector<int> vertices_to_labels = q.vertices_to_labels;
		vertices_to_labels.push_back(extended_label);

		/* Check DFS-code of candidate pattern in the map of discovered patterns */
		DFSCode *candidate_ = candidate_->GlobalMin(candidate_graph, vertices_to_labels);
		string candidate_DFSCode = candidate_->dfscode_str;
		cerr << "Candidate DFS Code: " << candidate_DFSCode << "\n";
		if (dictionary.count(candidate_DFSCode))
		{ /* Skip if candidate already in explored dictionary. 
			(NOTE: candidate frequency has not yet been computed) */
			global_dict_hits++;
			delete candidate_;
			priority_queue<subgraph_pattern, vector<subgraph_pattern>, ComparePattern> temp_leaf;
			int same_graph_in_leaf = 0;
			while (!(leaf.empty()))
			{
				subgraph_pattern qu = leaf.top();
				if (qu.strdfscode == candidate_DFSCode)
				{
					qu.subgraphs.insert(q.id);
					qu.subgraphs.insert(q.subgraphs.begin(), q.subgraphs.end());
					same_graph_in_leaf++;
				}
				leaf.pop();
				temp_leaf.push(qu);
			}
			assert(same_graph_in_leaf == 1);
			leaf = temp_leaf;
			continue;
		}

		/* Since pattern is absent in dictionary, construct replica and get support: */
		cerr << "GET_REPLICA call:\n";
		replica_s *replica = get_replica_brute(q, vert_id, extended_label, edge_label);
		cerr << "---------------------------\n";
		assert(q.replica->mappings.size() == replica->mappings.size() - 1);
		pair<int, int> sup_ind = support_index(replica);
		int &support = sup_ind.first;
		int &min_index = sup_ind.second;
		if (support >= min_sup)
		{
			/*Store constructed candidate as a subgraph_pattern: */
			vector<int> rm_path = candidate_->right_path();
			delete candidate_;
			subgraph_pattern candidate(global_num_patterns, candidate_graph, vertices_to_labels, replica, support, min_index, candidate_DFSCode, rm_path);
			dictionary.insert(candidate.strdfscode);
			candidate.subgraphs.insert(q.id);
			candidate.subgraphs.insert(q.subgraphs.begin(), q.subgraphs.end());
			candidate.subgraphs.insert(freq_vert_patternind.begin(), freq_vert_patternind.end());
			patternSizeToCount2[candidate.graph.size()] += 1;
			write_pattern(candidate.graph, candidate.vertices_to_labels, frequent_txt, global_num_patterns, candidate_DFSCode);
			leaf.push(candidate);
			global_num_patterns++;
		}
		else
		{
			cerr<<"deleting replica\n";
			delete replica;
		}
	}
	delete q.replica;
}

void possible_extensions_(subgraph_pattern &q, vector<LL> &all_extensions)
{
	//check if need to write a dfs code for a vertex
	cerr << "possible extensions called:\n";
	replica_s *replica = q.replica;
	vector<int> &rightmost_path = q.rightmost_path;
	unordered_set<LL> extensions;
	/*
	cerr << " printing Q:\n";
	for (int i = 0; i < q.graph.size(); i++)
	{
		cerr << i << ":";
		for (int j = 0; j < q.graph[i].size(); j++)
		{
			cerr << q.graph[i][j].first << ",";
		}
		cerr << "\n";
	}
	cerr << "DFS code is:" << q.strdfscode << "\n";
	cerr << "support is:" << q.support << "\n";
	sleep(pause_);
	*/
	for (int i = 0; i < rightmost_path.size(); i++)
	{
		unordered_set<int> &mappings = replica->mappings[rightmost_path[i]];
		vector<struct edge *> *dfslist = q.dfsordercons(rightmost_path[i]);
		for (const auto &vert_id : mappings)
		{
			vector<pair<int, int>> &adj_list = inp_graph[vert_id];
			if (adj_list.size() == q.graph[rightmost_path[i]].size())
				continue;
			assert(adj_list.size() > q.graph[rightmost_path[i]].size());
			for (int j = 0; j < adj_list.size(); j++)
			{
				pair<int, int> &edge = adj_list[j];
				LL extension_index = (LL)vertices_to_labels[edge.first] * (int)(pow(10, 7) + 0.5) + edge.second * 1000 + rightmost_path[i];
				if (extensions.count(extension_index))
					continue;
				bool ans = instance_checkerDFS(q, dfslist, rightmost_path[i], vert_id, edge.first);
				if (ans)
				{
					assert(extension_index >= 0); //overflow-check
					extensions.insert(extension_index);
				}
			}
		}
	}
	copy(extensions.begin(), extensions.end(), back_inserter(all_extensions));
}

bool instance_checkerDFS(subgraph_pattern &q, vector<struct edge *> *dfslist, int &extending_index, const int &extending_index_graphID, int &conflict_ID)
{
	unordered_map<int, int> instance = {{extending_index, extending_index_graphID}}; // key: pattern-ID, value: replica-ID
	set<int> instance_set = {extending_index_graphID};								 //simply to store all values of instance map elements for access
	int list_index = 1;																 //start DFS-wise instance enumeration from dfslist[1] (dfslist[0] is <extending_index, -1>)
	return DFS_instance(q, list_index, dfslist, instance, instance_set, conflict_ID);
}

bool DFS_instance(subgraph_pattern &q, int list_index, vector<struct edge *> *dfslist, unordered_map<int, int> &instance, set<int> &instance_set, int &conflict_ID)
{
	/*returns true if at least one instance found*/
	/*base-case:*/
	if (list_index == dfslist->size() || dfslist->size() == 1)
	{
		//dfslist has been completely traversed successfully
		assert(instance.size() == dfslist->size());
		return true;
	}

	/*recursion:*/
	for (const auto &child : q.replica->adj_list[instance[dfslist->at(list_index)->parent]])
	{
		//"child" iterates over all valid mappings of children of the replica vertex chosen as the mappping of its parent in the (possible) instance under consideration
		if (q.replica->vmaplist[child.first].count(dfslist->at(list_index)->child) && (child.second == dfslist->at(list_index)->edge_label))
		{
			//valid mapping of child found
			if (child.first == conflict_ID || instance_set.count(child.first))
				continue;
			/*child cannot be same as conflict_ID (ID of the attempted extended_index)
			Also, child cannot be a vertex already assigned as a mapping to some other pattern vertex of this instance*/
			instance[dfslist->at(list_index)->child] = child.first;
			instance_set.insert(child.first);
			if (DFS_instance(q, list_index + 1, dfslist, instance, instance_set, conflict_ID))
			{
				return true;
			}
			else
			{
				instance.erase(dfslist->at(list_index)->child);
				instance_set.erase(child.first);
			}
		}
	}
	return false;
}

pair<int, int> support_index(replica_s *replica)
{
	// returns the pair<MNI support, corresponding index> for the pattern from its replica
	int min_mapping_size = INT_MAX;
	int index = -1;
	for (int i = 0; i < replica->mappings.size(); i++)
	{
		if (replica->mappings[i].size() < min_mapping_size)
		{
			min_mapping_size = replica->mappings[i].size();
			index = i;
		}
	}
	pair<int, int> temp = make_pair(min_mapping_size, index);
	assert(index != -1);
	return temp;
}

/* Set-A: CSM-approximate replica-extension: */
replica_s *get_replica_brute(subgraph_pattern &q, int extending_index, int extended_label, int edge_label)
{
	assert(extending_index >= 0 && extended_label >= 0 && edge_label >= 0);
	cerr << "Replica call: extending index| extended label| edge label::" << extending_index << " " << extended_label << " " << edge_label << "\n";
	cerr << "top-k size:" << top_k.size() << "/" << k << "\n";
	cerr << "Parent support: " << q.support << "\n";
	for (int index = 0; index < q.graph.size(); index++)
	{
		cerr << q.vertices_to_labels[index] << " :";
		for (auto &neighbour : q.graph[index])
		{
			cerr << q.vertices_to_labels[neighbour.first] << "(" << neighbour.second << ")"
				 << ",";
		}
		cerr << "\n";
	}

	/*All structures: */
	//A. for parent replica:
	replica_s *parent_replica = q.replica;
	vector<vector<pair<int, int>>> *parent_graph = &q.graph;
	vector<unordered_set<int>> *parent_mappings = &q.replica->mappings;
	vector<struct edge *> *dfslist = q.dfsordercons(extending_index);
	vector<vector<int>> &parent_to_child = parent_replica->parent_to_child;
	parent_to_child.clear();
	parent_to_child.resize(parent_graph->size());
	for (int i = 1; i < dfslist->size(); i++)
	{
		parent_to_child[dfslist->at(i)->parent].push_back(dfslist->at(i)->child);
	}
	vector<int> &leaves = parent_replica->leaf_vertexID;
	leaves.clear();
	for (int i = 0; i < parent_to_child.size(); i++)
	{
		if (!parent_to_child[i].size())
		{
			leaves.push_back(i);
		}
	}
	parent_replica->enumerated_mappings.resize(parent_graph->size());
	
	//B. for candidate replica ("new replica"):
	replica_s *new_replica = new replica_s(parent_replica->adj_list, parent_graph->size());
	vector<unordered_set<int>> &new_replica_mappings = new_replica->mappings;
	new_replica_mappings.resize(q.graph.size() + 1); //NOTE: new_replica_mappings size is +1 original pattern size for storing extension mappings.
	
	/* Iterating extending index mappings: */
	int ext_siz = 0;
	bool found = false; //valid-extension
	for (const auto &mapped_extending : (*parent_mappings)[extending_index])
	{
		if (parent_graph->size() >= 3)
			cerr << "\r" << ++ext_siz << "/" << (*parent_mappings)[extending_index].size();
		pattern_sizes.push_back((*parent_graph).size());
		pfile << (*parent_graph).size() << "\n";
		found = false;
		for (const auto &input_nbr : inp_graph[mapped_extending])
		{
			if ((vertices_to_labels[input_nbr.first] == extended_label) && (input_nbr.second == edge_label))
			{
				found = false; //Is set to true if at least one valid instance found with the chosen "extension edge mapping"
				/*Enumerate all instances (using DFS) with the chosen "extension edge mapping":*/
				instance_allDFS(q, dfslist, extending_index, mapped_extending, input_nbr.first, new_replica_mappings, found); 
				if (found) 
				{	/*At least one valid instance found with the chosen "extension edge mapping"*/
					new_replica_mappings[parent_graph->size()].insert(input_nbr.first);
					new_replica->adj_list[mapped_extending].push_back(input_nbr);
					new_replica->adj_list[input_nbr.first].push_back(make_pair(mapped_extending, edge_label));					
				}
			}
		}
	}

	/*Assign new_replica_mappings to current new_replica mappings and update vmaplist*/
	for (int i = 0; i < new_replica_mappings.size(); i++)
	{
		for (const auto &set_elem : new_replica_mappings[i])
		{
			new_replica->vmaplist[set_elem].insert(i);
		}
	}

	/*Reset q: first visit map, potential child map, confirmed child map*/
	parent_replica->first_visit.erase(parent_replica->first_visit.begin(), parent_replica->first_visit.end());
	parent_replica->potential_child_mapping.erase(parent_replica->potential_child_mapping.begin(), parent_replica->potential_child_mapping.end());
	parent_replica->confirmed_child_mapping.erase(parent_replica->confirmed_child_mapping.begin(), parent_replica->confirmed_child_mapping.end());
	parent_replica->enumerated_mappings.clear();
	return new_replica;
}

void instance_allDFS(subgraph_pattern &q, vector<struct edge *> *dfslist, int &extending_index, const int &extending_index_vertexID, const int &conflict_ID, vector<unordered_set<int>> &new_mappings, bool &found)
{
	unordered_map<int, int> instance; // key: pattern_vtxID, value: replica_vtxID (=input_graph_vtxID)
	instance[extending_index] = extending_index_vertexID;
	set<int> instance_set = {extending_index_vertexID};
	int list_index = 1; //start DFS-wise instance enumeration from dfslist[1] (dfslist[0] is <extending_index, -1>)
	DFS_instance_all(q, list_index, dfslist, instance, instance_set, conflict_ID, new_mappings, found);
}

int DFS_instance_all(subgraph_pattern &q, int list_index, vector<struct edge *> *dfslist, unordered_map<int, int> &instance, set<int> &instance_set, const int &conflict_ID, vector<unordered_set<int>> &new_mappings, bool &global_found)
{
	/*returns true iff at least one instance found*/

	bool local_found = false; //"true" iff including current child_mapping results in at least one valid instance
	/* base case*/
	if (list_index == dfslist->size())
	{
		//dfslist has been traversed till the end
		assert(instance.size() == dfslist->size());
		for (auto &instance_mapping : instance)
		{
			new_mappings[instance_mapping.first].insert(instance_mapping.second);
		}
		for (const auto &leaf_vertex : q.replica->leaf_vertexID)
		{
			q.replica->enumerated_mappings[leaf_vertex].insert(instance[leaf_vertex]);
		}
		global_found = true;
		return local_found = true;
	}
	
	/* recursion */
	/*If first-visit then store set of all possible child-mapping vertices, else iterate over existing set of remaining potential child mappings*/
	int &parent_patternID = dfslist->at(list_index)->parent;
	int &parent_vertexID = instance[parent_patternID];
	int &child_patternID = dfslist->at(list_index)->child;
	int &edge_label = dfslist->at(list_index)->edge_label;
	unordered_set<int> &potential_set = q.replica->potential_child_mapping[to_string(parent_vertexID) + ' ' + to_string(child_patternID)];
	unordered_set<int> &confirmed_set = q.replica->confirmed_child_mapping[to_string(parent_vertexID) + ' ' + to_string(child_patternID)];

	if (!q.replica->first_visit[to_string(parent_patternID) + ' ' + to_string(child_patternID)].count(parent_vertexID))
	{ 
		//this is the first-visit at pattern_vertexID for this edge-type
		for (const auto &edge : q.replica->adj_list[parent_vertexID])
		{
			if (q.replica->vmaplist[edge.first].count(child_patternID) && edge.second == edge_label) //&& edge.first!=conflict_ID)
			{
				potential_set.insert(edge.first);
				/*store the set of potentially-valid child mappings*/
			}
		}
		q.replica->first_visit[to_string(parent_patternID) + ' ' + to_string(child_patternID)].insert(parent_vertexID);
	}

	/*iterate over (remaining) potential child mappings and try to enumerate instances*/
	unordered_set<int> add_to_confirmed;
	for (auto &child_mapping : potential_set)
	{
		if (instance_set.count(child_mapping) || child_mapping == conflict_ID)
			continue;
		instance[child_patternID] = child_mapping;
		instance_set.insert(child_mapping);
		bool instanceFound = DFS_instance_all(q, list_index + 1, dfslist, instance, instance_set, conflict_ID, new_mappings, global_found);
		if (instanceFound)
			local_found = true;
		/*Check if child_mapping is "completely-enumerated"*/
		if (instanceFound && q.replica->enumerated_mappings[child_patternID].count(child_mapping))
		{
			add_to_confirmed.insert(child_mapping);
			//to be deleted from potential_child_mapping and move to confirmed_child_mapping
		}
		instance.erase(child_patternID);
		instance_set.erase(child_mapping);
	}

	if (!local_found)
	{
		/*in this case (i.e. not a single valid instance found using potential set), 
		iterate over confirmed child mappings and try to enumerate (one) instance*/
		for (auto &child_mapping : confirmed_set)
		{
			if (instance_set.count(child_mapping) || child_mapping == conflict_ID)
				continue;
			instance[child_patternID] = child_mapping;
			instance_set.insert(child_mapping);
			bool instanceFound = DFS_instance_all(q, list_index + 1, dfslist, instance, instance_set, conflict_ID, new_mappings, global_found);
			instance.erase(child_patternID);
			instance_set.erase(child_mapping);
			if (instanceFound)
			{
				local_found = true;
				break;
			}
		}
	}
	
	/*iterate over add_to_confirmed set and update potential_child_mapping (deletion) and confirmed_child_mapping (insertion) */
	confirmed_set.insert(add_to_confirmed.begin(), add_to_confirmed.end());
	for (auto &child_mapping : add_to_confirmed)
	{
		potential_set.erase(child_mapping);
	}
	add_to_confirmed.clear();
	
	/*if ALL left-over children sets are empty, parent-vertex is also completely enumerated FOR THE CORRESPONDING MAPPING */
	bool ToBeMarkedEnumerated = true;
	for (auto &child_pattern : q.replica->parent_to_child[parent_patternID])
	{
		string key = to_string(parent_vertexID) + ' ' + to_string(child_pattern);
		if (!q.replica->potential_child_mapping.count(key) || q.replica->potential_child_mapping[key].size())
		{
			ToBeMarkedEnumerated = false;
			break;
		}
	}
	if (ToBeMarkedEnumerated)
	{
		q.replica->enumerated_mappings[parent_patternID].insert(parent_vertexID);
	}
	return local_found;
	/* 
	for (const auto &child : q.replica->adj_list[instance[dfslist[list_index].second]])
	{
		//"child" iterates over all valid mappings of children of the replica vertex chosen as the mappping of its parent in the (possible) instance under consideration
		if (q.replica->vmaplist[child].count(dfslist[list_index].first))
		{
			//possibly valid mapping of child found
			if (child == conflict_ID || instance_set.count(child))
				continue;
			instance[dfslist[list_index].first] = child;
			instance_set.insert(child);
			DFS_instance_all(q, list_index + 1, dfslist, instance, instance_set, conflict_ID, new_mappings, found);
			instance.erase(dfslist[list_index].first);
			instance_set.erase(child);
		}
	}
	*/
}

/* Set-B: CSM-complete replica-extension: */
replica_s *get_replica_brute_orig(subgraph_pattern &q, int extending_index, int extended_label, int edge_label)
{
	/* Extension of form A->E. "A#": A marked for retention. 
    Some metadata prints: */
	std::clock_t start;
	double duration;
	start = std::clock();
	cerr << "BRUTE replica call: extending index| extended label| edge label::" << extending_index << " " << extended_label << " " << edge_label << "\n";
	cerr << "Q:" << q.support << "\n";
	cerr << "top-k size:" << top_k.size() << "/" << k << "\n";
	int index = -1;
	for (auto &pattern_vtx : q.graph)
	{
		index++;
		cerr << q.vertices_to_labels[index] << " :";
		for (auto &neighbour : pattern_vtx)
		{
			cerr << q.vertices_to_labels[neighbour.first] << "(" << neighbour.second << ")"
				 << ",";
		}
		cerr << "\n";
	}

	/* All structures: */
	replica_s *parent_replica = q.replica;
	vector<vector<pair<int, int>>> *Q_graph = &q.graph;
	vector<unordered_set<int>> *Q_mappings = &q.replica->mappings;
	bool found = false; //valid-extension
	replica_s *new_replica = new replica_s(parent_replica->adj_list, Q_graph->size());
	vector<unordered_set<int>> new_replica_mappings(q.graph.size() + 1); //will replace current new_replica->mappings at the end
	//NOTE: new_replica_mappings size is +1 original pattern size. For extension mappings.
	cout << "new_replica initialisation size: " << new_replica->adj_list.size() << " new_replica_mappings size: " << new_replica_mappings.size() << "\n";
	cout << "new_replica_mappings.size prev:" << new_replica_mappings.size() << "\n";
	/* Iterating extending index mappings: */
	int ext_siz = 0;
	/*initialise first_visit members*/
	vector<struct edge *> *dfslist = q.dfsordercons(extending_index);
	cerr << "dfslist BRUTE_ORIG:\n";
	for (auto &elem : *dfslist)
	{
		cerr << elem->parent << "-->" << elem->child << "=" << elem->edge_label << "\n";
	}

	vector<vector<int>> parent_to_child(q.graph.size(), vector<int>());
	for (int i = 1; i < dfslist->size(); i++)
	{
		unordered_set<int> temp;
		q.replica->first_visit[to_string(dfslist->at(i)->parent) + '$' + to_string(dfslist->at(i)->child)] = temp;
		parent_to_child[dfslist->at(i)->parent].push_back(dfslist->at(i)->child);
	}
	q.replica->parent_to_child = parent_to_child;
	vector<int> leaves;
	for (int i = 0; i < parent_to_child.size(); i++)
	{
		if (!parent_to_child[i].size())
		{
			leaves.push_back(i);
		}
	}
	q.replica->leaf_vertexID = leaves;
	vector<unordered_set<int>> enumerated_mappings(q.graph.size(), unordered_set<int>());
	q.replica->enumerated_mappings = enumerated_mappings;
	for (const auto &mapped_extending : (*Q_mappings)[extending_index])
	{
		if (q.graph.size() >= 3)
			cerr << "\r" << ext_siz++ << "/" << (*Q_mappings)[extending_index].size();
		pattern_sizes.push_back((*Q_graph).size());
		pfile << (*Q_graph).size() << "\n";
		found = false;
		for (auto &input_nbr : inp_graph[mapped_extending])
		{
			if ((vertices_to_labels[input_nbr.first] == extended_label) && (input_nbr.second == edge_label))
			{
				found = false;
				instance_allDFS_orig(q, dfslist, extending_index, mapped_extending, input_nbr.first, new_replica_mappings, found);
				if (found)
				{
					new_replica_mappings[new_replica_mappings.size() - 1].insert(input_nbr.first);
					new_replica->adj_list[mapped_extending].push_back(input_nbr);
					if (new_replica->adj_list.count(input_nbr.first))
						new_replica->adj_list[input_nbr.first].push_back(make_pair(mapped_extending, edge_label));
					else
					{
						vector<pair<int, int>> adjlist_newvtx = {make_pair(mapped_extending, edge_label)};
						// set<int> adjlist_newvtx = {mapped_extending};
						new_replica->adj_list[input_nbr.first] = adjlist_newvtx;
					}
				}
			}
		}
	}
	cout << "new_replica_mappings.size after:" << new_replica_mappings.size() << "\n";
	/* Assign new_replica_mappings to current new_replica mappings and update vmaplist*/
	new_replica->mappings = new_replica_mappings;
	for (auto &map_iterator : new_replica->vmaplist)
	{
		set<int> empty;
		map_iterator.second = empty;
	}
	for (int i = 0; i < new_replica_mappings.size(); i++)
	{
		cout << "new_replica_mappings[i].size:" << new_replica_mappings[i].size() << "\n";
		for (const auto &set_elem : new_replica_mappings[i])
		{
			new_replica->vmaplist[set_elem].insert(i);
		}
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	cerr << "time: " << duration << '\n';
	return new_replica;
}

void instance_allDFS_orig(subgraph_pattern &q, vector<struct edge *> *dfslist, int &extending_index, int extending_index_vertexID, int &conflict_ID, vector<unordered_set<int>> &new_mappings, bool &found)
{
	unordered_map<int, int> instance; // <patternvtxID, replicavtxID>
	instance[extending_index] = extending_index_vertexID;
	set<int> instance_set = {extending_index_vertexID};
	int list_index = 1; //start DFS-wise instance enumeration from dfslist[1] (dfslist[0] is <extending_index, -1>)
	DFS_instance_all_orig(q, list_index, dfslist, instance, instance_set, conflict_ID, new_mappings, found);
}

void DFS_instance_all_orig(subgraph_pattern &q, int list_index, vector<struct edge *> *dfslist, unordered_map<int, int> &instance, set<int> &instance_set, int &conflict_ID, vector<unordered_set<int>> &new_mappings, bool &found)
{
	/* base case*/
	if (list_index == dfslist->size())
	{
		//dfslist has been traversed till the end
		assert(instance.size() == dfslist->size());
		for (auto &instance_mapping : instance)
		{
			if (instance.size() >= 3 && (instance[2] == 4045 || instance[2] == 4101))
				cerr << instance_mapping.first << "::" << instance_mapping.second << "\n";
			new_mappings[instance_mapping.first].insert(instance_mapping.second);
		}
		found = true;
		return;
	}
	/* recursion */
	/*If not first-visit then store set of all possible child-mapping vertices, else iterate over existing set of remaining potential child mappings*/
	int parent_patternID = dfslist->at(list_index)->parent;
	int parent_vertexID = instance[parent_patternID];
	int child_patternID = dfslist->at(list_index)->child;
	int edge_label = dfslist->at(list_index)->edge_label;

	for (const auto &child : q.replica->adj_list[instance[dfslist->at(list_index)->parent]])
	{
		/*"child" iterates over all valid mappings of children of the replica vertex, 
		chosen as the mappping of its parent in the (possible) instance under consideration*/
		if (q.replica->vmaplist[child.first].count(child_patternID) && (child.second == edge_label))
		{
			//possibly valid mapping of child found
			if (child.first == conflict_ID || instance_set.count(child.first))
				continue;
			instance[child_patternID] = child.first;
			instance_set.insert(child.first);
			DFS_instance_all_orig(q, list_index + 1, dfslist, instance, instance_set, conflict_ID, new_mappings, found);
			instance.erase(child_patternID);
			instance_set.erase(child.first);
		}
	}
}

void print_pattern(vector<vector<pair<int, int>>> &graph, vector<int> &vertices_to_labels)
{
	cout << "Pattern:\n";
	for (int i = 0; i < graph.size(); i++)
	{
		cout << i << ": ";
		for (auto &nbr : graph[i])
			cout << nbr.first << "(" << nbr.second << ") ";
		cout << "\n";
	}
	for (int i = 0; i < graph.size(); i++)
	{
		cerr << vertices_to_labels[i] << ": ";
		for (auto &nbr : graph[i])
			cerr << vertices_to_labels[nbr.first] << " ";
		cerr << "\n";
	}
}

void write_pattern(vector<vector<pair<int, int>>> &graph, vector<int> &vertices_to_labels, ofstream &out_, int id = -1, string dfscode_ = "_")
{
	if (id != -1)
		out_ << "Pattern ID: " << id << " ";
	if (dfscode_ != "_")
		out_ << "DFScode: " << dfscode_ << "\n";
	out_ << "Adjacency list:\n";
	for (int i = 0; i < graph.size(); i++)
	{
		out_ << i << ": ";
		for (auto &nbr : graph[i])
			out_ << nbr.first << "(" << nbr.second << ") ";
		out_ << "\n";
	}
	out_ << "Labels:\n";
	for (int i = 0; i < graph.size(); i++)
	{
		out_ << vertices_to_labels[i] << ": ";
		for (auto &nbr : graph[i])
			out_ << vertices_to_labels[nbr.first] << " ";
		out_ << "\n";
	}
	out_ << "-----------------------------------------\n";
}