#include "csm_GS.h"

using namespace std;

void freq_vertices()
{
	frequent_vertices.resize(inp_graph.size(), false); //Initialize all vertices as NOT frequent
	cerr << "Label_to_vertices size: " << label_to_vertices.size() << "\n";
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

			/*Constructing and storing all instances for the frequent pattern:*/
			vector<vector<int>> instances;
			
			for (int j=0; j<label_to_vertices[i].size(); j++)
			{
				frequent_vertices[label_to_vertices[i][j]]=true;
				vector<int> instance;
				instance.push_back(label_to_vertices[i][j]);
				instances.push_back(instance);
			}
			instances_s *inst = new instances_s(instances);

			/*Constructing and storing the replica for the frequent pattern:*/
			// vector<unordered_set<int>> mappings(1);
			// unordered_set<int> &maps = mappings[0];
			// unordered_map<int, set<edge_tuple>> replica_adj_list;
			// unordered_map<int, set<int>> vmaplist;
			// for (int j = 0; j < label_to_vertices[i].size(); j++) //iterating over all vertices of the frequent label[i]
			// {
			// 	frequent_vertices[label_to_vertices[i][j]] = true;
			// 	vmaplist[label_to_vertices[i][j]] = {0}; //pattern-vtxID is trivially zero
			// 	replica_adj_list[label_to_vertices[i][j]] = {};
			// 	maps.insert(label_to_vertices[i][j]); // "maps": stores all vtxID's of freq label[i]
			// }
			// replica_s *replica_ = new replica_s(mappings, replica_adj_list, vmaplist); //"replica_": replica structure for frequent vertex

			/*Constructing the frequent pattern graph and the subgraph_pattern data structure:*/
			vector<vector<pair<int, int>>> graph_(1);
			vector<int> vertices_to_labels = {i};
			int support = label_to_vertices[i].size();
			int min_index = 0;
			string dfscode_ = to_string(i); //DFS-code of vertices is the vertex label
			vector<int> rightmost_path = {0};
			subgraph_pattern pattern_(id, graph_, vertices_to_labels, inst, support, min_index, dfscode_, rightmost_path);
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
	//NOTE: only deletion from input graph adjacency list; no vertex deletion to maintain original vertex-ID's
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


void search()
{
	int count = 0;
	while (!leaf.empty())
	{
		clock_t start = clock();
		//cout<<"search counter: "<<++search_counter<<", q: ";
		// cerr << "Ceasing condition not met. Search loop counter: " << ++count
		// 	 << "\nLeaf Size: " << leaf.size() << ", Top-k size: " << top_k.size() << "/" << k << "\n";
		subgraph_pattern q = leaf.top();
		cout<<q.strdfscode<<", post-op: "<<endl;
		cout<<"-------------------"<<endl;
		leaf.pop();
		patternSizeToCount[q.graph.size()] += 1;
		operated.push_back(q);
		// operate(q);
		//double project_time = (double)(clock() - project_start_time) / CLOCKS_PER_SEC;
		//cout<<project_time<<", post-ex: ";
		
		extend(q);
		cerr << "Post extend() Leaf size:" << leaf.size() << "\n";
		//project_time = (double)(clock() - project_start_time) / CLOCKS_PER_SEC;
		//cout<<project_time<<"s, top-k-size: "<<top_k.size()<<"\n";
	}
}



void extend(subgraph_pattern &q)
{
	cerr << "extend() call:\n";
	vector<vector<pair<int, int>>> &Q_graph = q.graph;
	vector<int> &Q_vertices_to_labels = q.vertices_to_labels;
	unordered_set<LL> extensions;
	possible_extensions_(q, extensions);
	cout<<" extensions-size: "<<extensions.size()<<" "; 
	int extensions_size = extensions.size();
	cerr << "Candidate extensions size: " << extensions_size << "\n";
	pair<int, int> sup_ind_org = support_index(q.instances);
	int &support_org = sup_ind_org.first;
	int &min_index_org = sup_ind_org.second;
	// for (int i = 0; i < extensions.size(); i++)
	// {
	// 	cerr << "Candidate extn #" << i + 1 << ": " << extensions[i] << "\n";
	// }

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

	//for (int i = 0; i < extensions_size; i++)
	for (const auto& extension: extensions)
	{
		//cerr << "extension: " << i + 1 << "/" << extensions_size << ": ";
		LL extension_ = extension;
		int vert_id = extension_ % 1000;
		extension_ /= 1000;
		int edge_label = extension_ % 10000;
		extension_ /= 10000;
		int extended_label = extension_;
		cerr << " vert_id: " << vert_id << " edge_l: " << edge_label << " extended_l: " << extended_label << "\n";

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
			cout<<same_graph_in_leaf<<endl;
			assert(same_graph_in_leaf <= 1);
			leaf = temp_leaf;
			continue;
		}

		instances_s *instances;
		instances = get_instances(q, vert_id, extended_label, edge_label);



		cerr << "---------------------------\n";
		// assert(q.replica->mappings.size() == replica->mappings.size() - 1);
		pair<int, int> sup_ind = support_index(instances);
		int &support = sup_ind.first;
		int &min_index = sup_ind.second;
		if (support >= min_sup)
		{
			/*Store constructed candidate as a subgraph_pattern: */
			vector<int> rm_path = candidate_->right_path();
			delete candidate_;
			subgraph_pattern candidate(global_num_patterns, candidate_graph, vertices_to_labels, instances, support, min_index, candidate_DFSCode, rm_path);
			dictionary.insert(candidate.strdfscode);
			candidate.subgraphs.insert(q.id);
			candidate.subgraphs.insert(q.subgraphs.begin(), q.subgraphs.end());
			candidate.subgraphs.insert(freq_vert_patternind.begin(), freq_vert_patternind.end());
			// write_pattern(candidate.graph, candidate.vertices_to_labels, output_txt, global_num_patterns, candidate_DFSCode);
			output_txt << "(rm-path: ";
			for (int j = 0; j < rm_path.size(); j++)
				output_txt << rm_path[j] << " ";
			output_txt << ")\n";
			leaf.push(candidate);
			global_num_patterns++;
		}
		else
		{
			cerr << "deleting replica\n";
			delete instances;
		}
	}
	delete q.instances;
}

void possible_extensions_(subgraph_pattern &q, unordered_set<LL> &all_extensions)
{
	cerr << "possible extensions called:\n";
	instances_s *instances_curr = q.instances;
	vector<int> &rightmost_path = q.rightmost_path;
	vector<vector<int>> instances = instances_curr->instances;
	int instance_number = instances.size();
	//unordered_set<LL> extensions;
	for (int i = 0; i < rightmost_path.size(); i++)
	{
		
		for(int k=0;k<instance_number;k++){
			vector<struct edge *> *dfslist = q.dfsordercons(rightmost_path[i]);
			int vert_id = instances[k][rightmost_path[i]];
			vector<pair<int, int>> &adj_list = inp_graph[vert_id];
			if (adj_list.size() == q.graph[rightmost_path[i]].size())
				continue;

			assert(adj_list.size() > q.graph[rightmost_path[i]].size());
			vector<int> curr_inst = instances[k];
			for (int j = 0; j < adj_list.size(); j++)
			{
				pair<int, int> &edge = adj_list[j];
				LL extension_index = (LL)vertices_to_labels[edge.first] * (int)(pow(10, 7) + 0.5) + edge.second * 1000 + rightmost_path[i];
				if (all_extensions.count(extension_index))
					continue;
				if(find(curr_inst.begin(), curr_inst.end(), edge.first)==curr_inst.end()){
					assert(extension_index >= 0); //overflow-check
					//cout<<extension_index<<"\n";
					all_extensions.insert(extension_index);
				}
			}
		}
		
	}
	//copy(extensions.begin(), extensions.end(), back_inserter(all_extensions));
}



pair<int, int> support_index(instances_s *instances){
	unordered_map<int, unordered_set<int>> ma;
	int si = instances->instances[0].size();
	unordered_set<int> tem;
	for(int i=0;i<si;i++){
		ma[i] = tem;
	}
	for(int i=0;i<instances->instances.size();i++){
		for(int j=0;j<si;j++){
			ma[j].insert(instances->instances[i][j]);
		}
	}
	int min_mapping_size = INT_MAX;
	int index = -1;
	for(int i=0;i<si;i++){
		if(ma[i].size()<min_mapping_size){
			min_mapping_size = ma[i].size();
			index = i;
		}
	}
	pair<int, int> temp = make_pair(min_mapping_size, index);
	return temp;
}
 
instances_s *get_instances(subgraph_pattern &q, int extending_index, int extended_label, int edge_label){
	vector<vector<int>> instances;
	vector<vector<int>> old_instances = q.instances->instances;
	int si = old_instances.size();
	for(int i=0;i<si;i++){
		vector<int> inst_curr = old_instances[i];
		int vert_id = inst_curr[extending_index];
		vector<pair<int, int>> adj_list = inp_graph[vert_id];
		int siz_adj = adj_list.size();
		for(int j=0;j<siz_adj;j++){
			if(vertices_to_labels[adj_list[j].first]==extended_label && adj_list[j].second==edge_label && find(inst_curr.begin(), inst_curr.end(), adj_list[j].first)==inst_curr.end()){
				vector<int> instance_n = inst_curr;
				instance_n.push_back(adj_list[j].first);
				instances.push_back(instance_n);
			}
		}
	}
	instances_s *inst_ret = new instances_s(instances);
	return inst_ret;

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
	cout << "Labels:\n";
	for (int i = 0; i < graph.size(); i++)
	{
		cout << vertices_to_labels[i] << ": ";
		for (auto &nbr : graph[i])
			cout << vertices_to_labels[nbr.first] << " ";
		cout << "\n";
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
	out_ << "---------------------------------------------------------------------\n";
}
