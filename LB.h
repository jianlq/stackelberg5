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

#endif