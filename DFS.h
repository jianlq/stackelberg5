#ifndef DFS_H
#define DFS_H
#include"CGraph.h"
#include"EE.h"

void CGraph::KSP(int s, int t, unsigned int k)
{
	//FILE *out = fopen("outputFile//ksp.txt", "a");
	listSolu.clear();
	listTemp.clear();
	Status *beginning = new Status(s,0);
	beginning->passby.push_back(s);
	listTemp.push_back(beginning);
	//fprintf(out, "ksp while 111   "); 
	while(listSolu.size() < k && listTemp.size() != 0)
	{
		listTemp.sort(pStatusComp);///在更新序列里进行排序
		int id = (**listTemp.begin()).ID;///拿出节点序列的第一个点的ID
		if(id == t)
			listSolu.push_back(*listTemp.begin());///如果这个点就是终点，那么把这个状态丢到解集合里头去。
		double d = (**listTemp.begin()).d;///拿出存储的d值
		list<int> passby = (**listTemp.begin()).passby;
		vector<CEdge*>::iterator it;
		for(it = mapVID_listEdge[id].begin(); it != mapVID_listEdge[id].end(); it++)///遍历该点的出度边
		{
			int head = (**it).getHead();///拿出头结点ID
			if(find(passby.begin(),passby.end(),head) == passby.end())
			{
				double weight = (**it).getWeight() + d;///拿出该边的重量和d相加得到路径长度
				list<int> templist = passby;
				templist.push_back(head);
				Status *newstatus = new Status(head,weight,*listTemp.begin(),templist);///形成某个状态
				listTemp.push_back(newstatus);///把状态丢进更新序列里
			}
		}
		listTemp.pop_front();
	}
	//fprintf(out, "ksp while 222   "); 
	listPath.clear(); ////////每个OD对之间算的路径集合  必须清空
	list<Status*>::iterator Sit;
	for(Sit = listSolu.begin(); Sit != listSolu.end(); Sit++){
		listPath.push_back(new CPath(beginning, *Sit, mapVID_Vertex, mapVID_listEdge));
	}
	//fprintf(out, "ksp while 333\n");  
	//fclose(out);
}


//计算每个OD对之间的路径
////  计算每个 OD对之间的路径   
bool CGraph::GAinit(vector<demand> &req){
	
	////KSP
	//reqlistPath.clear();
	//for(unsigned int i=0;i<req.size();i++){  
	//	KSP(req[i].org,req[i].des,K); //计算得到一个OD对的listPath
	//	reqlistPath.push_back(listPath);
	//	if(listPath.size() < K )
	//		return false;
	//}
	//return true;

	////// DFS
	reqlistPath.clear();
	for(unsigned int i=0;i<req.size();i++){
		int j = 0;
		vector<vector<CEdge*>> reqpath;
		while(j < K){
			SetUNVISITED();
			pathver.clear();
			myDFS(req[i].org,req[i].des);
			if(DFSflag){
				j++;
				vector<CEdge*> tmp;
				int p = 0;
				while( p < (pathver.size()-1) ){
					int s = pathver[p],t = pathver[p+1];
					for(unsigned int  a = 0; a <Link.size(); a++){
						if(Link[a]->tail == s && Link[a]->head == t ){
							tmp.push_back(Link[a]);
							break;
						}
					}
					p++;
				}
				reqpath.push_back(tmp);
			}
		}
		if( reqpath.size() < K)
			return false;
		reqlistPath.push_back(reqpath);		
	}
	return true;
}

void CGraph::myDFS(int s,int t){
		if (!visit[s]){
		   pathver.push_back(s);
		   visit[s] = true;
		}
		if( s == t ){   // 找到了终点
			DFSflag = 1;
			return ; 
		}
		vector<int> unvit;
		int ss ,ret=0; 
		for(unsigned  int i = 0; i < adjL[s].size();i++){
			if(visit[adjL[s][i]->head])
				continue;	
			else
				unvit.push_back(adjL[s][i]->head); //没有访问的顶点编号		
		}		
		if(unvit.size()== 0) {  // 没有未被访问的邻节点
			pathver.pop_back();
			if(pathver.size()>0) { 
				ss = pathver[pathver.size()-1];
				ret = 1;  //表示有回溯点
			}
		}
		else{
			int j = rand() % unvit.size();
			ss = unvit[j]; 
			ret = 1;  //表示有继续深入点
		}
		if(ret)
			myDFS(ss,t);
		else
			return;
 }


//load balance
void heuristicLB(CGraph *G,vector<demand>&req,int ornum,double &mlu,double &delay){
	G->clearOcc();
	double util = 0;
	int block = 0,success = 0;
	for(unsigned int i = 0; i < req.size(); i++){
		double ret = G->dijkstra(i,req[i].org, req[i].des, req[i].flow,1,0,1); 
		if(ret >= INF)
			block++;
		else{
			success++;
			util = max(mlu,ret);
		}
	}
	mlu = INF,delay = INF;
	if(success == req.size()){	
		mlu = util;

		//mlu
		mlu = 0;
		for(int i=0;i<G->m;i++){  
			mlu = max(mlu,G->Link[i]->use/G->Link[i]->capacity);
		}

		//delay
		delay = 0;
		for(int i=0;i<G->m;i++){  
			G->Link[i]->latency = linearCal(G->Link[i]->use,G->Link[i]->capacity);
		}
		for(int d=0;d < ornum;d++){
			double cur = 0;
			for(unsigned int ij=0;ij<G->reqPathID[d].size();ij++)
				cur +=  G->Link[G->reqPathID[d][ij]]->latency;
			if(LINEAR)
				delay += cur;
			else
				delay += cur*req[d].flow;
		}

	}	
}

// delay
void heuristicOR(CGraph *G,vector<demand>&req,int ornum,double &delay,double &mlu){
	G->clearOcc();
	int block = 0,success = 0;
	for(unsigned int i = 0; i < req.size(); i++){
		double ret = G->dijkstra(i,req[i].org, req[i].des, req[i].flow,0,1,1); 
		if(ret+1e-5 >= INF)
			block++;
		else{
			success++;		
		}
	}
	mlu = delay = INF;
	if(success == req.size()){	
		mlu = 0,delay = 0;
		//mlu
		mlu = 0;
		for(int i=0;i<G->m;i++){  
			mlu = max(mlu,G->Link[i]->use/G->Link[i]->capacity);
		}

		//delay
		delay = 0;
		for(int i = 0; i < G->m; i++){  
			G->Link[i]->latency = linearCal(G->Link[i]->use,G->Link[i]->capacity);
		}
		for(int d = 0;d < ornum;d++){
			double cur = 0;
			for(unsigned int ij = 0;ij < G->reqPathID[d].size(); ij++)
				cur +=  G->Link[G->reqPathID[d][ij]]->latency;
			if(LINEAR)
				delay += cur;
			else
				delay += cur*req[d].flow;
		}
		
		cout << "HOR   "<<mlu<<"\t"<<delay<<endl;
	}
}


#endif