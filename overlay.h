#ifndef OVERLAY
#define OVERLAY
#include"Common.h"
#include <ilcplex/ilocplex.h>
#include"CGraph.h"

#define epsilon (double)0.001
extern double alpha;
extern double beta;


class overlay
{
	private:
	int No;
	int n;
	int m;
	int reqnum;
	double cost_value;
	vector<int> mapping; //点映射到下层
	vector<vector<int>> adj_mtx;
	vector<double> delay;
public:
	vector<demand> req;
	vector<vector<double>> flow_mark;//用来标记每个流怎么分
	vector<demand> edge_flow;

	public:
	overlay()
	{
		No = NULL;
		n = 0;
		m = 0;
		reqnum = 0;
		cost_value = INF;
		mapping.clear();
		adj_mtx.clear();
		delay.clear();
		req.clear();
		flow_mark.clear();
		edge_flow.clear();
	}
	overlay(int number,string filename)
	{
		reqnum = 0;
		No = number;
		cost_value = INF;
		mapping.clear();
		adj_mtx.clear();
		delay.clear();
		req.clear();
		flow_mark.clear();
		edge_flow.clear();
		
		ifstream fin(filename);
		fin>>n;
		for(int i = 0;i < n;i++)
		{
			int temp;
			fin>>temp;
			mapping.push_back(temp);
		}
		fin>>m;

		vector<int> adj_row;
		for(int i = 0;i < m;i++)
			adj_row.push_back(-1);
		for(int i = 0;i < m;i++)
			adj_mtx.push_back(adj_row);


		for(int i = 0;i < m;i++)
		{
			delay.push_back(0.1);
			int tail,head;
			fin>>tail>>head;
			adj_mtx[tail][head] = i;//// tail -----> head
			demand dem(mapping[tail],mapping[head],epsilon);
			edge_flow.push_back(dem);
		}
	}

	void show()
	{
		cout<<"show overlay"<<No<<":"<<endl;
		cout<<n<<"  "<<m<<endl;
		for(int i = 0;i < n;i++)
			cout<<mapping[i]<<" ";
		cout<<endl;
		for(int i = 0;i < n;i++)
		{
			for(int j = 0;j < n;j++)
				cout<<adj_mtx[i][j]<<"  ";
			cout<<endl;
		}
		for(int i = 0;i < m;i++)
			cout<<edge_flow[i].org<<"  "<<edge_flow[i].des<<" "<<edge_flow[i].flow<<endl;
	}
	vector<demand>& getTrafficMatrix()
	{
		return edge_flow;
	}
	int getNo()
	{
		return No;
	}

	int getEdgeNum(){
		return m;
	}

	void get_OD_infor(string filename)
	{
		ifstream fin(filename);
		fin>>reqnum;
		req.clear();
		for(int i = 0;i < reqnum;i++)
		{
			int src,dst,flow;
			fin>>src>>dst>>flow;
			demand dem(src,dst,flow);
			req.push_back(dem);
		}
	}

	void showOD(){
		for(int i = 0;i < reqnum;i++)
		{
			cout<<req[i].org<<"  "<<req[i].des<<"  "<<req[i].flow<<endl;
		}
	}

	double LP()
	{
		IloEnv env;
		IloModel mod(env);
		IloCplex solver(mod);

		IloArray<IloIntVarArray> x(env, reqnum);
		for(int i = 0; i < reqnum; i++)
			x[i] = IloIntVarArray(env, m, 0, 1);
			
		//点的流量平衡
		for(int i = 0;i < reqnum;i++){
			for(int j = 0;j < n;j++)
			{
				IloExpr constraint(env);
				for(int k = 0;k < n;k++)
				{
					if(adj_mtx[j][k] != -1)
						constraint += x[i][adj_mtx[j][k]];
					if(adj_mtx[k][j] != -1)
						constraint -= x[i][adj_mtx[k][j]];
				}

				if(j == req[i].org)
					mod.add(constraint == 1);  ////////?????
				else if(j == req[i].des)
					mod.add(constraint == -1);
				else
					mod.add(constraint == 0);
			}
		}

		//cost 
		IloExpr cost(env);
		for(int i = 0;i < m;i++)
		{
			IloExpr constraint(env);
			for(int j = 0;j < reqnum;j++)
				constraint += x[j][i];                //// Linear fitting
				//constraint += x[j][i] * req[j].flow; //// x1/(C-x)   1/C-x :delay[i]
			cost += delay[i] * constraint; 
		}
		mod.add(IloMinimize(env, cost));

		solver.setOut(env.getNullStream());
		double obj = INF;
		if(solver.solve()){
			obj = cost_value = solver.getObjValue();
			cout <<"overlay "<<No<<" cost "<<obj<<endl;

			if(flow_mark.size() == 0)
			{
				for(int i = 0;i < reqnum;i++)
				{
					vector<double> mark_row;
					for(int j = 0;j < m;j++)
					{
						mark_row.push_back(solver.getValue(x[i][j]));
					}
					flow_mark.push_back(mark_row);
				}
			}
			else
			{
				for(int i = 0;i < reqnum;i++)
				{

					for(int j = 0;j < m;j++)
					{
						flow_mark[i][j] = (double)(1-beta)*flow_mark[i][j] + (double)beta * solver.getValue(x[i][j]);
					}
				}
			}

			for(int i = 0;i < m;i++)
			{
				double f = 0;
				for(int j = 0;j < reqnum;j++)
					f += flow_mark[j][i] * req[j].flow;
				f += epsilon;
				edge_flow[i].flow = f;
			}
		}
		else{
			cout << "overlay "<< No <<" unfeasible"<<endl;
		}
		env.end();
		return obj;

	}
	void getDelay(vector<double> delay_from_underlying)
	{
		delay = delay_from_underlying;
	}
	
	void updateDelay(vector<double> delay_from_underlying)
	{
		for(int i = 0;i < m;i++)
			delay[i] = (1-alpha)*delay[i] + alpha*delay_from_underlying[i];
	}


	void printMark()
	{
		cout<<"<<<<<<<<<<<<<mark<<<<<<<<<<<<<<<"<<endl;
		for(int i = 0;i < reqnum;i++)
		{
			for(int j = 0; j< m;j++)
				cout<<flow_mark[i][j]<<"  ";
			cout<<endl;
		}
	}
	void logValue(string filename)
	{
		ofstream fout(filename,ios::app);
		fout<<cost_value<<endl;
	}
	double getCost()
	{
		return cost_value;
	}
	
	void printDelay()
	{
		for(int i = 0;i < delay.size();i++)
			cout<<delay[i]<<endl;
	}
	void clearDelay(){
		for(int i = 0;i < delay.size();i++)
			delay[i] = 0.1;
	}
};

#endif
