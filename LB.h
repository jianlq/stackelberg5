#ifndef LB_H
#define LB_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

// 4+x^2
double LBdictor(CGraph *G,vector<demand> & req){
	G->clearOcc();
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	//优化目标
	IloNumVar z(env, 0, 1);	
	model.add(IloMinimize(env, z));

	// 对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];

			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		for(int d = 0; d <  num; d++)
			load += req[d].flow*x[d][i];
		model.add( load <= G->Link[i]->capacity );  
		model.add( load <= z*G->Link[i]->capacity);	
	}	

	EEsolver.setOut(env.getNullStream());
	double obj = INF;
	if(EEsolver.solve()){
		obj = EEsolver.getObjValue(); //mlu
		//delay
		for(int i=0;i<G->m;i++){  
			double loadc = 0;
			for(int d=0;d < num;d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->latency = linearCal(loadc,G->Link[i]->capacity);
		}
		double del = 0;
		for(int d=0;d < num; d++){  //ornum
			for(int i=0;i<G->m;i++)
				del +=  EEsolver.getValue(x[d][i])*G->Link[i]->latency;
		}
		G->delay = del;
		cout << "LB\t利用率 "<< obj <<"\t延时 "<<del<<endl;
	}
	else{
		cout << "LB unfeasible"<<endl;
	}
	for(int i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
	return obj;
}

double ORdictor(CGraph *G,vector<demand>& req,int start,int end){    
	G->clearOcc();
	IloEnv env;
	IloModel model(env);
	IloCplex ORsolver(model);
	cout <<start<<" start end "<<end<<endl;
	int totalnum = req.size();
	IloArray<IloNumVarArray> x(env, totalnum); 
	for(int d = 0; d < totalnum; d++)
		x[d] = IloNumVarArray(env, G->m, 0, 1); // num * G->m 的二维矩阵

	IloNumVarArray cost(env,G->m,0.0 ,INF); 
	IloExpr res(env);
	for(int i=0;i<G->m;i++)
	{
		IloNumExpr load(env);
		for(int d = start;d < end; d++)
			load += x[d][i]*req[d].flow;
		double c = G->Link[i]->capacity;
		model.add(cost[i] >= load);
		model.add(cost[i] >= 3*load-2*c/3);
		model.add(cost[i] >= 10*load-16*c/3);
		model.add(cost[i] >= 70*load-178*c/3);
		res += cost[i];
	}
	model.add(IloMinimize(env, res));

	// 流量守恒约束
	for(int d = 0; d < totalnum; d++)
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}

		//带宽约束
		for(int i = 0; i < G->m; i++){
			IloExpr constraint(env);
			for(int d = 0; d <  totalnum; d++)
				constraint += req[d].flow*x[d][i];
			model.add(constraint <= G->Link[i]->capacity);  
		}

		ORsolver.setOut(env.getNullStream());
		double obj = INF;

		if(ORsolver.solve()){
			obj = ORsolver.getObjValue();
			//	cout << obj <<endl;
			// mlu
			double util = 0;
			for(int i = 0; i < G->m; i++){  
				double load = 0;
				for(int d = 0; d < totalnum; d++)
					load += ORsolver.getValue(x[d][i])*req[d].flow;
				util = max(util,load/G->Link[i]->capacity);
				G->Link[i]->use = load;
			}
			G->mlu = util;

			// delay
			double del = 0;
			for(int i = 0; i < G->m; i++){  
				double cur = linearCal(G->Link[i]->use,G->Link[i]->capacity);
				for(int d = start; d < end; d++)
					del += ORsolver.getValue(x[d][i])*cur;
			}
			obj = del;
			cout << "OR\t利用率 "<<util<< "\t延时 "<<obj<<endl;
		}
		else
			env.out() << "OR unfeasible" << endl;

		for(int i = 0; i < req.size(); i++)
			x[i].end();
		x.end();
		env.end();
		return obj;
}

void TE(CGraph *G,vector<demand> &background,vector<vector<demand>> &overlay)
{
	int x_num = 0;
	for(int i = 0;i < overlay_num;i++)
		x_num += overlay[i].size(); /////x_num表示总的overlay的需求个数

	IloEnv env;
	IloModel model(env);
	IloCplex solver(model);

	//变量
	IloArray<IloNumVarArray> x(env, x_num); ////x代表overlay的分配
	for(int i = 0; i < x_num; i++)
		x[i] = IloNumVarArray(env, G->m, 0, 1);

	int background_num = background.size();
	IloArray<IloNumVarArray> y(env, background_num); ////y代表背景流的分配
	for(int i = 0; i < background_num; i++)
		y[i] = IloNumVarArray(env, G->m, 0, 1);


	IloNumVar z(env,0,1);
	model.add(IloMinimize(env, z));

	//边的容量限制
	for(int i = 0;i < G->m; i++)
	{
		IloExpr constraint(env);
		for(int j = 0;j < background_num;j++)
			constraint += y[j][i] * background[j].flow;

		int count = 0;
		for(int j = 0;j < overlay_num;j++)
		{
			for(int k = 0; k < overlay[j].size(); k++)
			{
				constraint += x[count][i] * overlay[j][k].flow;
				count++;
			}
		}
		model.add(constraint <= z * G->Link[i]->capacity);
	}
	//点的背景流量平衡
	for(int d = 0; d < background_num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];
			// 出 - 入
			if(i == background[d].org)
				model.add(constraint == 1);
			else if(i == background[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	//点的overlay流量平衡
	int count = 0;
	for(int i = 0;i < overlay_num;i++)
	{
		for(int j = 0;j < overlay[i].size();j++)
		{
			for(int k = 0;k < G->n; k++) // nodes
			{
				IloExpr constraint(env);
				for(unsigned int l = 0; l < G->adjL[k].size(); l++) // 出度边 link
					constraint += x[count][G->adjL[k][l]->id];
				for(unsigned int l = 0; l < G->adjRL[k].size(); l++) // 入度边 link
					constraint -= x[count][G->adjRL[k][l]->id];
				if(k == overlay[i][j].org)
					model.add(constraint == 1);
				else if(k == overlay[i][j].des)
					model.add(constraint == -1);
				else
					model.add(constraint == 0);
			}
			count++;
		}
	}

	solver.setOut(env.getNullStream());
	if(solver.solve()){
		for(int i = 0;i < background_num;i++)
			for(int j = 0;j < G->m;j++)
				G->background_mark[i][j] = solver.getValue(y[i][j]);

		cout<<"obj value:"<<solver.getObjValue()<<endl;

		count = 0;
		for(int i = 0;i< overlay_num;i++)
		{
			for(int j = 0;j < overlay[i].size();j++)
			{	
				for(int k = 0; k < G->m ; k++)
					G->overlay_mark[i][j][k] = solver.getValue(x[count][k]);
				count++;
			}
		}

		//////////更新delay
		for(int i = 0;i < G->m; i++)
		{
			double flow = 0;

			for(int j = 0;j < background_num;j++)
				flow += G->background_mark[j][i] *background[j].flow;

			for(int j = 0;j < overlay_num;j++)
			{
				for(int k = 0; k < overlay[j].size(); k++)
				{
					flow += G->overlay_mark[j][k][i] * overlay[j][k].flow;
				}
			}
			//G->Link[i]->latency =1.0/(G->Link[i]->capacity- flow); // 1/(C-x)
			G->Link[i]->latency = linearCal(flow, G->Link[i]->capacity); //linear fitting
		}
		//////////////更新 to_overlay
		for(int i = 0;i < overlay_num;i++)
		{
			for(int j = 0;j < overlay[i].size();j++)
			{
				float sum = 0;
				for(int k = 0; k < G->m; k++)
					sum += G->overlay_mark[i][j][k] * G->Link[k]->latency;
				G->to_overlay[i][j] = sum ;
			}

		}
	}
	else{
		cout<<"TE unfeasible"<<endl;
	}
	env.end();

}

#endif