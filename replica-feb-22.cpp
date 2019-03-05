#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iterator>
#include "csm.h"

using namespace std;

unordered_map<int, int> patternSizeToCount;
unordered_map<int, int> patternSizeToCount2;

unsigned int pause_ = PAUSE_TIME;

vector<int> pattern_sizes;
unordered_set<string> dictionary;
vector<set<int> *> corvm;
vector<set<int>> corP;
int min_index_change_ctr = 0;
int global_counter = 0;
ofstream pfile;
ofstream frequent_txt;

vector<int> vertices_to_labels;
int global_num_patterns = 0;
vector<bool> frequent_vertices;
int k = 20;
int min_sup;
int control;
int hop;
vector<int> freq_vert_patternind;
int global_dict_hits = 0;

vector<vector<int>> label_to_vertices; //label to vertices
vector<vector<pair<int, int>>> inp_graph;

priority_queue<subgraph_pattern, vector<subgraph_pattern>, ComparePattern> leaf;
vector<subgraph_pattern> operated;
priority_queue<correlated_patterns, vector<correlated_patterns>, CompareCorrelation> top_k;

int main(int argc, char **argv)
{
	clock_t project_start_time = clock();
	frequent_txt.open("frequent_patterns.txt", ios::out);
	pfile.open("pfile.txt", ios::out);
	string input_file = argv[1];
	cout << "Running top-K CSM for graph: " << input_file << "\n";
	min_sup = atoi(argv[2]);
	cout << "min_sup: " << min_sup << "\n";
	hop = atoi(argv[5]);
	cout << "hop: " << hop << "\n";
	cout << "K: " << k << "\n";
	ofstream filw(argv[3]);
	cout << "writing to file: " << argv[3] << "\n";
	control = atoi(argv[4]);
	cout << "control: " << control << "\n";
	cout << "***********start***********\n";

	/*Parsing input file to obtain vertex count and the largest label:*/
	clock_t measure_time_start = clock();
	ifstream inFile;
	inFile.open(input_file);
	if (!inFile)
	{
		cout << "Unable to open file:" << input_file;
		exit(1);
	}
	int vertices = 0, maxLabel = -1, line_num = 1, vertex_label, vertexnum, vertexnum2, edgelabel, edgeind;
	string x;
	while (getline(inFile, x))
	{
		istringstream iss(x);
		vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
		if (tokens[0][0] == 'v')
		{
			vertices++;
			stringstream(tokens[2]) >> vertex_label;
			if (vertex_label > maxLabel)
				maxLabel = vertex_label;
		}
		else if (tokens[0][0] == 'e')
		{
			break;
		}
		else
		{
			cout << "Incorrect input format. Line:" << line_num << ": " << x << "\n";
			exit(1);
		}
		line_num++;
	}
	inFile.close();
	label_to_vertices.resize(maxLabel + 1);
	inp_graph.resize(vertices);
	corP.resize(vertices);
	cerr << "First-parse time:" << (double)(clock() - measure_time_start) / CLOCKS_PER_SEC << "s \n";

	/*Parsing input file to load the input graph:*/
	inFile.open(input_file);
	if (!inFile)
	{
		cout << "Unable to open file:" << input_file;
		exit(1);
	}
	while (getline(inFile, x))
	{
		istringstream iss(x);
		vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
		if (tokens[0][0] == 'v')
		{
			stringstream(tokens[2]) >> vertex_label;
			vertices_to_labels.push_back(vertex_label);
			stringstream(tokens[1]) >> vertexnum;
			label_to_vertices[vertex_label].push_back(vertexnum);
		}
		else if (tokens[0][0] == 'e')
		{
			stringstream(tokens[1]) >> vertexnum;
			stringstream(tokens[2]) >> vertexnum2;
			stringstream(tokens[3]) >> edgelabel;
			assert(edgelabel<10000); //Edge-label should necessarily be less the 1.E+04
			inp_graph[vertexnum].push_back(make_pair(vertexnum2, edgelabel));
			inp_graph[vertexnum2].push_back(make_pair(vertexnum, edgelabel));
		}
	}
	inFile.close();
	cerr << "Total parse time: " << (double)(clock() - measure_time_start) / CLOCKS_PER_SEC << "s \n\n";

	/*Obtain all frequent labels and corresponding vertices:*/
	measure_time_start = clock();
	freq_vertices();
	cerr << "F-1 time: " << (double)(clock() - measure_time_start) / CLOCKS_PER_SEC << "s \n\n";

	/*Modify adjacency lists based on freq_vertices() results:*/
	measure_time_start = clock();
	modify_adjlist();
	cerr << "Adjacency -list modification time: " << (double)(clock() - measure_time_start) / CLOCKS_PER_SEC << "s \n\n";
	cerr << "UNIX Pause_\n";
	//sleep(pause_);

	/*CorV computation follows infrequent vertices-based adj-list modification:*/
	measure_time_start = clock();
	compute_corv(inp_graph, corvm, hop);
	assert(corvm.size() == inp_graph.size());
	cout << "CorV computation time: " << (double)(clock() - measure_time_start) / CLOCKS_PER_SEC << "s\n\n";

	/*search() controls CSM program execution hereafter:*/
	measure_time_start = clock();
	search();
	cout << "search() execution time: " << (double)(clock() - measure_time_start) / CLOCKS_PER_SEC << "\n\n";

	cout << "Computation Finished. Number of operated patterns: "<<operated.size()<<"\n";
	
	cerr<<"Remaining Leaf supports: ";
	while (!leaf.empty())
	{
		const subgraph_pattern &q = leaf.top();
		cerr<<q.support<<" ";
		leaf.pop();
	}
	cerr<<"\n";

	//print all operated patterns:
	int op_id = 0;
	for (auto &pattern : operated)
	{
		cerr << "Operated pattern #" << op_id++ << "/" << operated.size() <<" "<<pattern.strdfscode<<" \n";
		for (int i = 0; i < pattern.graph.size(); i++)
		{
			cerr << pattern.vertices_to_labels[i] << ":";
			for (int j = 0; j < pattern.graph[i].size(); j++)
			{
				pair<int, int> &temp = pattern.graph[i][j];
				cerr << pattern.vertices_to_labels[temp.first] << "(" << temp.second << ") ";
			}
			cerr << "\n";
		}
	}

	cout << top_k.size() << " " << top_k.top().corr_value << endl;

	while (!top_k.empty())
	{
		for (int i = 0; i<operated.size(); i++)
		{
			subgraph_pattern &qu = operated[i];
			if (qu.graph.size()>1)
			{
				frequent_txt<<i<<":"<<qu.id<<", ";
			}
		}
		frequent_txt<<"\n";
		const correlated_patterns &corr_pair = top_k.top();
		subgraph_pattern &q1 = operated[corr_pair.operated_id_1];
		subgraph_pattern &q2 = operated[corr_pair.operated_id_2];
		cerr << "Correlated Pattern IDs:" << q1.id << ", " << q2.id << " correlation value:" << corr_pair.corr_value << "\n";
		cerr << "Pattern 1:\n";
		for (int i = 0; i < q1.graph.size(); i++)
		{
			cerr << q1.vertices_to_labels[i] << ":";
			for (int j = 0; j < q1.graph[i].size(); j++)
			{
				pair<int, int> temp = q1.graph[i][j];
				cerr << q1.vertices_to_labels[temp.first] << "(" << temp.second << ") ";
			}
			cerr << "\n";
		}
		cerr<<q1.strdfscode<<"\n";
		// 
		// for (int i = 0; i < q.vertices_to_labels.size(); i++)
		// {
		// 	cout << q.vertices_to_labels[i] << endl;
		// }
		// 
		cerr << "Pattern 2:\n";
		for (int i = 0; i < q2.graph.size(); i++)
		{
			cerr << q2.vertices_to_labels[i] << ":";
			for (int j = 0; j < q2.graph[i].size(); j++)
			{
				pair<int, int> temp = q2.graph[i][j];
				cerr << q2.vertices_to_labels[temp.first] << "(" << temp.second << ") ";
			}
			cerr << "\n";
		}
		cerr<<q2.strdfscode<<"\n";
		// 
		// for (int i = 0; i < q.vertices_to_labels.size(); i++)
		// {
		// 	cout << q.vertices_to_labels[i] << endl;
		// }
		// 
		top_k.pop();
	}
	
	filw.close();

	cerr << "min index change ctr:" << min_index_change_ctr << "\n";

	cerr << "global ctr:" << global_counter << "\n";
	pfile.close();

	cerr << "PatternSizeToCount:\n";
	for (auto &elem : patternSizeToCount)
	{
		cerr << elem.first << ":" << elem.second << "\n";
	}
	cerr << "p2sc2:\n";
	for (auto &elem : patternSizeToCount2)
	{
		cerr << elem.first << ":" << elem.second << "\n";
	}
	pfile.close();
	frequent_txt.close();
	cerr<<"Total dictionary hits: "<<global_dict_hits<<"\n";
	double project_time = (double)(clock() - project_start_time) / CLOCKS_PER_SEC;
	cout << "Total project runnning time: " << project_time << " seconds.\n";
}
