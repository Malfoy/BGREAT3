/*****************************************************************************
 *   Bgreat : De Bruijn Graph Read Mapping Tool
 *   Authors: Antoine Limasset
 *   Contact: antoine.limasset@gmail.com
 *   Source: https://github.com/Malfoy/BGREAT
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty ofF
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *consen
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/



#include <thread>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>



uint trimingBases(0);



//~ vector<uNumber> Aligner::alignReadGreedy(const string& read, bool& overlapFound, uint errors, bool& rc){
	//~ vector<pair<kmer,uint>> listOverlap(getNOverlap(read,tryNumber));
	//~ if(listOverlap.empty()){++noOverlapRead;return {};}
	//~ overlapFound=true;
	//~ vector<uNumber> pathBegin,pathEnd;
	//~ for(uint start(0); start<(uint)listOverlap.size(); ++start){
		//~ pathBegin={};
		//~ uint errorBegin(checkBeginGreedy(read,listOverlap[start],pathBegin,errors));
		//~ if(errorBegin<=errors){
			//~ pathEnd={};
			//~ uint errorsEnd(checkEndGreedy(read,listOverlap[start],pathEnd,errors-errorBegin));
			//~ if(errorsEnd+errorBegin<=errors){
				//~ reverse(pathBegin.begin(),pathBegin.end());
				//~ pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
				//~ return pathBegin;
			//~ }
		//~ }
	//~ }
	//~ return {};
//~ }


double Aligner::alignment_weight(const vector<int>& path){
	double res(0);
	double base(k-1);
	for(uint i(0);i<path.size();++i){
		res+=unitigs_weight[abs(path[i])]*(unitigs[abs(path[i])].size()-k+1);
		base+=(unitigs[abs(path[i])].size()-k+1);
	}
	return res/base;
}



//OKOK
alignment Aligner::alignReadGreedyAnchors(const string& read, int score_max, const pair<pair<uint,uint>,uint>& anchor){
	//~ cout<<"arga"<<endl;
	alignment al;
	string unitig("");
	bool returned(false);
	int score;
	int unitigNumber(anchor.first.first),positionUnitig(anchor.first.second),positionRead(anchor.second);
	alignment al_begin(al),al_end(al);
	if(unitigNumber>=0){
		unitig=(unitigs[unitigNumber]);
	}else{
		if(rcMode){
			unitig=(unitigsRC[-unitigNumber]);
		}else{
			unitig=reverseComplements(unitigs[-unitigNumber]);
		}
		positionUnitig=unitig.size()-positionUnitig-anchorSize;
		returned=true;
	}
	al.position_anchors_in_read=positionRead;
	al.position_anchors_in_unitig=positionUnitig;
	if(positionRead>positionUnitig){
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//CASE 1 : unitig included in read
			//~ cout<<"1:"<<endl;
			//TODO OPTIM START BY THE SHORTEST PIECE
			int local_score=(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,score_max));
			score=local_score;
			if(score>=score_min){
				score+=(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},al_begin,score_min));//TODO
				if(score>=score_min){
					//~ cout<<"1!"<<endl;
					score+=(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},al_end,score_max-score));//TODO
					//~ cout<<score<<endl;
					if(score>=score_max){//TODO
						//~ cout<<"1!!"<<endl;
						reverse(al_begin.path.begin(),al_begin.path.end());
						reverse(al_begin.scores.begin(),al_begin.scores.end());
						al.path=al_begin.path;
						al.scores=al_begin.scores;
						al.path.push_back(unitigNumber);
						al.unitig_number_anchored=al.path.size()-1;
						al.scores.push_back(local_score);
						al.path.insert(al.path.end(), al_end.path.begin(),al_end.path.end());
						al.scores.insert(al.scores.end(), al_end.scores.begin(),al_end.scores.end());
						al.score=score;
						return al;
					}
				}
			}
		}else{
			//CASE 2 : unitig overap read
			//~ cout<<"2:"<<endl;
			if(read.size()-positionRead+positionUnitig<k){
				//~ cout<<"nope"<<endl;
				al.score=-1000;
				return al;
			}
			int local_score=(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),score_max));
			score=local_score;
			if(score>=score_min){
				//~ cout<<"2.1:"<<endl;
				score+=(checkBeginGreedy(read,{str2num(unitig.substr(0,k-1)),positionRead-positionUnitig},al_begin,score_max-score));
				if(score>=score_max){
					//~ cout<<"2.2:"<<endl;
					reverse(al_begin.path.begin(),al_begin.path.end());
					reverse(al_begin.scores.begin(),al_begin.scores.end());
					al.path=al_begin.path;
					al.scores=al_begin.scores;
					al.path.push_back(unitigNumber);
					al.unitig_number_anchored=al.path.size()-1;
					al.scores.push_back(local_score);
					al.score=score;
					return al;
				}
			}
		}
	}else{
		if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//CASE 3 : read overlap unitig
			//~ cout<<"3:"<<endl;
			if(unitig.size()+positionRead-positionUnitig<k){
				//~ cout<<"nope"<<endl;
				al.score=-1000;
				return al;
			}
			//~ cout<<read<<endl;
			//~ cout<<read.substr(0,unitig.size()+positionRead-positionUnitig)<<endl;
			int local_score=(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),score_max));
			score=local_score;
			//~ cout<<unitig.size()+positionRead-positionUnitig<<" "<<score<<endl;
			if(score>=score_min){
				//~ cout<<"3.1:"<<endl;
				al.path.push_back(unitigNumber);
				al.scores.push_back(local_score);
				score+=(checkEndGreedy(read,{str2num(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},al_end,score_max-score));
				//~ cout<<score<<endl;
				if(score>=score_max){
					//~ cout<<"3.2:"<<endl;
					//~ cout<<score<<endl;
					al.unitig_number_anchored=0;
					al.path.insert(al.path.end(), al_end.path.begin(),al_end.path.end());
					al.scores.insert(al.scores.end(), al_end.scores.begin(),al_end.scores.end());
					al.score=score;
					return al;
				}
			}
		}else{
			//CASE 4 : read included in unitig
			//~ cout<<"4:"<<endl;
			score=(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,score_max));
			//~ cout<<unitig.substr(positionUnitig-positionRead,read.size())<<endl;
			//~ cout<<positionUnitig<<" "<<positionRead<<endl;
			//~ cout<<score<<endl;
			if(score>=score_max){
				//~ cout<<score<<" "<<score_max<<endl;
				al.unitig_number_anchored=0;
				al.path.push_back(unitigNumber);
				al.scores.push_back(score);
				al.score=score;
				return al;
			}else{
			}
		}
	}
	al.score=-1000;
	return al;
}


//~ vector<uNumber> Aligner::alignReadGreedyAnchorsstr(const string& read, uint errorMax, const pair<pair<uint,uint>,uint>& anchor, uint& errors){
	//~ vector<uNumber> pathBegin,pathEnd;
	//~ string unitig("");
	//~ bool returned(false);
	//~ int unitigNumber(anchor.first.first),positionUnitig(anchor.first.second),positionRead(anchor.second);
	//~ if(unitigNumber>=0){
		//~ unitig=(unitigs[unitigNumber]);
	//~ }else{
		//~ if(rcMode){
			//~ unitig=(unitigsRC[-unitigNumber]);
		//~ }else{
			//~ unitig=reverseComplements(unitigs[-unitigNumber]);
		//~ }
		//~ positionUnitig=unitig.size()-positionUnitig-anchorSize;
		//~ returned=true;
	//~ }

	//~ if(positionRead>=positionUnitig){
		//~ if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//~ errors=(missmatchNumber(read.substr(positionRead-positionUnitig,unitig.size()),unitig,errorMax));
			//~ if(errors<=errorMax){
				//~ pathBegin={};
				//~ errors+=(checkBeginGreedy(read,{(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				//~ if(errors<=errorMax){
					//~ pathEnd={(int)unitigNumber};
					//~ errors+=(checkEndGreedy(read,{(unitig.substr(unitig.size()-k+1,k-1)),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
					//~ if(errors<=errorMax){
						//~ reverse(pathBegin.begin(),pathBegin.end());
						//~ pathBegin.insert(pathBegin.end(), pathEnd.begin(),pathEnd.end());
						//~ return pathBegin;
					//~ }
				//~ }
			//~ }
		//~ }else{
			//~ errors=(missmatchNumber(read.substr(positionRead-positionUnitig),unitig.substr(0,read.size()-positionRead+positionUnitig),errorMax));
			//~ if(errors<=errorMax){
				//~ pathBegin={};
				//~ errors+=(checkBeginGreedy(read,{(unitig.substr(0,k-1)),positionRead-positionUnitig},pathBegin,errorMax-errors));
				//~ if(errors<=errorMax){
					//~ reverse(pathBegin.begin(),pathBegin.end());
					//~ pathBegin.push_back((int)unitigNumber);
					//~ return pathBegin;
				//~ }
			//~ }
		//~ }
	//~ }else{
		//~ if(read.size()-positionRead>=unitig.size()-positionUnitig){
			//~ errors=(missmatchNumber(unitig.substr(positionUnitig-positionRead),read.substr(0,unitig.size()+positionRead-positionUnitig),errorMax));
			//~ if(errors<=errorMax){
				//~ pathEnd={(int)positionUnitig-(int)positionRead,(int)unitigNumber};
				//~ errors+=(checkEndGreedy(read,{unitig.substr(unitig.size()-k+1,k-1),positionRead-positionUnitig+unitig.size()},pathEnd,errorMax-errors));
				//~ if(errors<=errorMax){
					//~ return pathEnd;
				//~ }
			//~ }
		//~ }else{
			//~ uint errors(missmatchNumber(unitig.substr(positionUnitig-positionRead,read.size()),read,errorMax));
			//~ if(errors<=errorMax){
				//~ return {(int)positionUnitig-(int)positionRead,(int)unitigNumber};
			//~ }
		//~ }
	//~ }
	//~ return {};
//~ }



//~ uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	//~ if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	//~ string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	//~ vector<pair<string,uNumber>> rangeUnitigs;
	//~ rangeUnitigs=(getEnd(overlap.first));
	//~ uint miniMiss(1000),miniMissIndice(9);
	//~ bool ended(false),equal(false);
	//~ int offset(0);
	//~ kmer nextOverlap(0);
	//~ for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ unitig=(rangeUnitigs[i].first);
		//~ //case the rest of the read is too small
		//~ if(unitig.size()-k+1>=readLeft.size()){
			//~ uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ path.push_back(unitig.size()-readLeft.size()-k+1);
				//~ return 0;
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ ended=true;
				//~ equal=false;
				//~ offset=unitig.size()-readLeft.size()-k+1;
			//~ }
		//~ }else{
			//~ //case the read is big enough we want to recover a true overlap
			//~ uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return mapOnLeftEndGreedy(read , path, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);

			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ kmer overlapNum(str2num(unitig.substr(0,k-1)));
				//~ if(miss<miniMiss){
					//~ ended=false;
					//~ equal=false;
					//~ miniMiss=miss;
					//~ miniMissIndice=i;
					//~ nextUnitig=unitig;
					//~ nextOverlap=overlapNum;
				//~ }
			//~ }
		//~ }
	//~ }

	//~ if (miniMiss<=errors and not equal){
		//~ path.push_back(rangeUnitigs[miniMissIndice].second);
		//~ if(ended){
			//~ path.push_back(offset);
			//~ return miniMiss;
		//~ }
		//~ return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	//~ }
	//~ if(equal){
		//~ return 1000;
	//~ }
	//~ return miniMiss;
//~ }



//~ uint Aligner::mapOnLeftEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors){
	//~ if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	//~ string unitig,readLeft(read.substr(0,overlap.second)),nextUnitig;
	//~ vector<pair<string,uNumber>> rangeUnitigs;
	//~ rangeUnitigs=(getEnd(overlap.first));
	//~ uint miniMiss(1000),miniMissIndice(9);
	//~ bool ended(false),equal(false);
	//~ int offset(0);
	//~ string nextOverlap;
	//~ for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ unitig=(rangeUnitigs[i].first);
		//~ //case the rest of the read is too small
		//~ if(unitig.size()-k+1>=readLeft.size()){
			//~ uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()), readLeft, errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ path.push_back(unitig.size()-readLeft.size()-k+1);
				//~ return 0;
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ ended=true;
				//~ equal=false;
				//~ offset=unitig.size()-readLeft.size()-k+1;
			//~ }
		//~ }else{
			//~ //case the read is big enough we want to recover a true overlap
			//~ uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()-(unitig.size()-k+1)), errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return mapOnLeftEndGreedy(read , path, {(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, errors);
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ ended=false;
				//~ equal=false;
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ nextUnitig=unitig;
				//~ nextOverlap=((unitig.substr(0,k-1)));
			//~ }
		//~ }
	//~ }
	//~ if (miniMiss<=errors and not equal){
		//~ path.push_back(rangeUnitigs[miniMissIndice].second);
		//~ if(ended){
			//~ path.push_back(offset);
			//~ return miniMiss;
		//~ }
		//~ return miniMiss+mapOnLeftEndGreedy(read , path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, errors-miniMiss);
	//~ }
	//~ if(equal){
		//~ return 1000;
	//~ }
	//~ return miniMiss;
//~ }



//~ uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<kmer, uint>& overlap , uint errors){
	//~ string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	//~ if(readLeft.size()<trimingBases){return 0;}
	//~ vector<pair<string,uNumber>> rangeUnitigs;
	//~ rangeUnitigs=getBegin(overlap.first);
	//~ uint miniMiss(1001), miniMissIndice(9);
	//~ bool ended(false),equal(false);
	//~ kmer nextOverlap(0);
	//~ for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ unitig=(rangeUnitigs[i].first);
		//~ //case the rest of the read is too small
		//~ if(unitig.size()-k+1>=readLeft.size()){
			//~ uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return 0;
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ ended=true;
				//~ equal=false;
			//~ }
		//~ }else{
			//~ //case the read is big enough we want to recover a true overlap
			//~ uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return (mapOnRightEndGreedy(read , path, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ ended=false;
				//~ equal=false;
				//~ kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ nextUnitig=unitig;
				//~ nextOverlap=overlapNum;
			//~ }
		//~ }
	//~ }
	//~ if(miniMiss<=errors and not equal){
		//~ path.push_back(rangeUnitigs[miniMissIndice].second);
		//~ if(ended){return miniMiss;}
		//~ return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	//~ }
	//~ if(equal){
		//~ return 1000;
	//~ }
	//~ return miniMiss;
//~ }



//~ uint Aligner::mapOnRightEndGreedy(const string &read, vector<uNumber>& path, const pair<string, uint>& overlap , uint errors){
	//~ string unitig,readLeft(read.substr(overlap.second)),nextUnitig;
	//~ if(readLeft.size()<trimingBases){return 0;}
	//~ vector<pair<string,uNumber>> rangeUnitigs;
	//~ rangeUnitigs=getBegin(overlap.first);
	//~ uint miniMiss(1000), miniMissIndice(9);
	//~ bool ended(false),equal(false);
	//~ string nextOverlap;
	//~ for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ unitig=(rangeUnitigs[i].first);
		//~ //case the rest of the read is too small
		//~ if(unitig.size()-k+1>=readLeft.size()){
			//~ uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()), readLeft, errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return 0;
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ ended=true;
				//~ equal=false;
			//~ }
		//~ }else{
			//~ //case the read is big enough we want to recover a true overlap
			//~ uint miss(missmatchNumber(unitig.substr(k-1), readLeft.substr(0,unitig.size()-k+1), errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return (mapOnRightEndGreedy(read , path, {(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)}, errors));
			//~ }else if(miss==miniMiss){
				//~ equal=true;
			//~ }else if(miss<miniMiss){
				//~ ended=false;
				//~ equal=false;
				//~ miniMiss=miss;
				//~ miniMissIndice=i;
				//~ nextUnitig=unitig;
				//~ nextOverlap=(unitig.substr(unitig.size()-k+1,k-1));
			//~ }
		//~ }
	//~ }
	//~ if(miniMiss<=errors and not equal){
		//~ path.push_back(rangeUnitigs[miniMissIndice].second);
		//~ if(ended){return miniMiss;}
		//~ return (miniMiss+mapOnRightEndGreedy(read , path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)}, errors-miniMiss));
	//~ }
	//~ if(equal){
		//~ return 1000;
	//~ }
	//~ return miniMiss;
//~ }


vector<uNumber> reverseVector(const vector<uNumber>& V){
	vector<uNumber> RC;
	for(uint i(0);i<V.size();++i){
		RC.push_back(-V[V.size()-i-1]);
	}
	return RC;
}


//~ void Aligner::path_extension_right(vector<uNumber>& path){
	//~ uNumber end(path[path.size()-1]);
	//~ string unitig;
	//~ unitig=getUnitig(end);
	//~ vector<pair<string,uNumber>> rangeUnitigs;
	//~ rangeUnitigs=(getBegin(str2num((unitig.substr(unitig.size()-k+1,k-1)))));
	//~ if(rangeUnitigs.size()==1){
		//~ path.push_back(rangeUnitigs[0].second);
	//~ }
//~ }



//~ void Aligner::path_extension(vector<uNumber>& path){
	//~ uint size(path.size());
	//~ path_extension_right(path);
	//~ path=reverseVector(path);
	//~ path_extension_right(path);
	//~ path=reverseVector(path);
	//~ if(path.size()!=size){
		//~ path_extension(path);
	//~ }
//~ }



//OKOK
int Aligner::checkBeginGreedy(const string& read,const pair<kmer, uint>& overlap, alignment& al, int score_max){
	//~ cout<<"CBG"<<endl;
	//~ cout<<score_max<<" "<<(int)overlap.second<<endl;
	if((int)overlap.second<score_max){return -1000;}
	//~ cout<<"CBG1"<<endl;
	string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getEnd(overlap.first));
	int best_score(-1000),indiceMinMiss(9);
	bool ended(false),equal(false);
	int offset(0);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ cout<<"CBG2"<<endl;

		unitig=rangeUnitigs[i].first;
		//~ cout<<unitig.size()<<endl;
		if(unitig.size()-k+1>=readLeft.size()){
			//~ cout<<"CBG2.1"<<endl;
			int score=(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft,score_min));
			if(score==readLeft.size()){
				//~ cout<<"CBG2.1.1"<<endl;
				al.path.push_back(rangeUnitigs[i].second);
				al.scores.push_back(score);
				return score;
			}else if(score==best_score){
				//~ cout<<"CBG2.1.2"<<endl;
				equal=true;
			}else if(score>best_score){
				//~ cout<<"CBG2.1.3"<<endl;
				best_score=score;
				indiceMinMiss=i;
				ended=true;
				equal=false;
				offset=unitig.size()-readLeft.size()-k+1;
			}
		}else{
			//~ cout<<"CBG2.2"<<endl;
			int score(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), score_min));
			if(score==unitig.size()-k+1){
				al.path.push_back(rangeUnitigs[i].second);
				al.scores.push_back(score);
				return score+checkBeginGreedy(read, {str2num(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)}, al,score_max-score);
			}else if(score==score_min){
				equal=true;
			}else if(score>best_score){
				kmer overlapNum(str2num(unitig.substr(0,k-1)));
				ended=false;
				equal=false;
				best_score=score;
				indiceMinMiss=i;
				nextUnitig=unitig;
				nextOverlap=overlapNum;
			}

		}
	}
	//~ cout<<"CBG3"<<endl;
	//~ if(best_score>0 and not equal){//TODO DECIDE WHAT TODO WHEN MULTIMAPPING
	if(best_score>score_min){
		al.path.push_back(rangeUnitigs[indiceMinMiss].second);
		al.scores.push_back(best_score);
		if(ended){
			//~ path.push_back(offset);
			return best_score;
		}
		return best_score+checkBeginGreedy(read,{nextOverlap,overlap.second-(nextUnitig.size()-k+1)}, al,score_max-best_score);
	}
	//~ if(equal){
		//~ return 1000;
	//~ }
	return best_score;
}


//TODO FOR STR
//~ uint Aligner::checkBeginGreedy(const string& read,const pair<string, uint>& overlap, vector<uNumber>& path, uint errors){
	//~ if(overlap.second<=trimingBases){path.push_back(0);return 0;}
	//~ string readLeft(read.substr(0,overlap.second)),unitig,nextUnitig;
	//~ vector<pair<string,uNumber>> rangeUnitigs;
	//~ rangeUnitigs=(getEnd(overlap.first));
	//~ uint minMiss(1001),indiceMinMiss(9);
	//~ bool ended(false),equal(false);
	//~ int offset(0);
	//~ string nextOverlap;
	//~ for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ unitig=rangeUnitigs[i].first;
		//~ if(unitig.size()-k+1>=readLeft.size()){
			//~ uint miss(missmatchNumber(unitig.substr(unitig.size()-readLeft.size()-k+1,readLeft.size()),readLeft, errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ path.push_back(unitig.size()-readLeft.size()-k+1);
				//~ return 0;
			//~ }else if(miss==minMiss){
				//~ equal=true;
			//~ }else if(miss<minMiss){
				//~ minMiss=miss;
				//~ indiceMinMiss=i;
				//~ ended=true;
				//~ equal=false;
				//~ offset=unitig.size()-readLeft.size()-k+1;
			//~ }
		//~ }else{
			//~ uint miss(missmatchNumber(unitig.substr(0,unitig.size()-k+1), readLeft.substr(readLeft.size()+k-1-unitig.size()), errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return mapOnLeftEndGreedy(read, path, {(unitig.substr(0,k-1)),overlap.second-(unitig.size()-k+1)},errors);
			//~ }else if(miss==minMiss){
				//~ equal=true;
			//~ }else if(miss<minMiss){
				//~ ended=false;
				//~ equal=false;
				//~ minMiss=miss;
				//~ indiceMinMiss=i;
				//~ nextUnitig=unitig;
				//~ nextOverlap=((unitig.substr(0,k-1)));
			//~ }

		//~ }
	//~ }
	//~ if(minMiss<=errors and not equal){
		//~ path.push_back(rangeUnitigs[indiceMinMiss].second);
		//~ if(ended){
			//~ path.push_back(offset);
			//~ return minMiss;
		//~ }
		//~ return minMiss+mapOnLeftEndGreedy(read, path, {nextOverlap,overlap.second-(nextUnitig.size()-k+1)},errors-minMiss);
	//~ }
	//~ if(equal){
		//~ return 1000;
	//~ }
	//~ return minMiss;
//~ }


//OKOK
int Aligner::checkEndGreedy(const string& read, const pair<kmer, uint>& overlap, alignment& al, int max_score){
	//~ cout<<"CEG"<<endl;

	string readLeft(read.substr(overlap.second)),unitig,nextUnitig;
	//~ cout<<(int)readLeft.size()<<endl;
	//~ cout<<max_score<<endl;
	//~ cout<<readLeft<<endl;
	if((int)readLeft.size()<max_score){return -1000;}
	//~ cout<<"CEG1"<<endl;
	vector<pair<string,uNumber>> rangeUnitigs;
	rangeUnitigs=(getBegin(overlap.first));
	int best_score(-1000),indiceMinMiss(9);
	bool ended(false),equal(false);
	kmer nextOverlap(0);
	for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ cout<<"CEG2"<<endl;
		unitig=rangeUnitigs[i].first;
		if(unitig.size()-k+1>=readLeft.size()){
			//~ cout<<"CEG1"<<endl;
			int score(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, score_min));
			//~ cout<<score<<endl;
			//~ cout<<readLeft.size()<<endl;
			if(score==readLeft.size()){
				al.path.push_back(rangeUnitigs[i].second);
				al.scores.push_back(score);
				return score;
			}else if(score==best_score){
				equal=true;
			}else if(score>best_score){
				best_score=score;
				indiceMinMiss=i;
				ended=true;
				equal=false;
			}
		}else{
			int score(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), score_min));
			//~ cout<<"CEG2"<<endl;
			//~ cout<<score<<endl;
			//~ cout<<unitig<<endl;
			//~ cout<<"noend"<<endl;
			if(score==unitig.size()-k+1){
				al.path.push_back(rangeUnitigs[i].second);
				al.scores.push_back(score);
				return score+checkEndGreedy(read, {str2num(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},  al,max_score-score);
			}else if(score==best_score){
				equal=true;
			}else if(score>best_score){
				kmer overlapNum(str2num(unitig.substr(unitig.size()-k+1,k-1)));
				best_score=score;
				indiceMinMiss=i;
				nextOverlap=overlapNum;
				nextUnitig=unitig;
				ended=false;
				equal=false;
			}
		}
	}
	//~ cout<<"CEG4"<<endl;
	if(best_score>score_min){
		//~ cout<<"CEG5"<<endl;
		al.path.push_back(rangeUnitigs[indiceMinMiss].second);
		al.scores.push_back(best_score);
		if(ended){return best_score;}
		return best_score+checkEndGreedy(read, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},al, max_score-best_score);
	}
	//~ if(equal){
		//~ return 1000;
	//~ }
	return best_score;
}


//TODO FOR STR
//~ uint Aligner::checkEndGreedy(const string& read, const pair<string, uint>& overlap, vector<uNumber>& path, uint errors){
	//~ string readLeft(read.substr(overlap.second)),unitig,nextUnitig;
	//~ if(readLeft.size()<=trimingBases){return 0;}
	//~ vector<pair<string,uNumber>> rangeUnitigs(getBegin(overlap.first));
	//~ uint minMiss(1000),indiceMinMiss(9);
	//~ bool ended(false),equal(false);
	//~ string nextOverlap;
	//~ for(uint i(0); i<rangeUnitigs.size(); ++i){
		//~ unitig=rangeUnitigs[i].first;
		//~ if(unitig.size()-k+1>=readLeft.size()){
			//~ uint miss(missmatchNumber(unitig.substr(k-1,readLeft.size()),readLeft, errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return 0;
			//~ }
			//~ else if(miss==minMiss){
				//~ equal=true;
			//~ }else if(miss<minMiss){
				//~ minMiss=miss;
				//~ indiceMinMiss=i;
				//~ ended=true;
				//~ equal=false;
			//~ }
		//~ }else{
			//~ uint miss(missmatchNumber(unitig.substr(k-1),readLeft.substr(0,unitig.size()-k+1), errors));
			//~ if(miss==0){
				//~ path.push_back(rangeUnitigs[i].second);
				//~ return mapOnRightEndGreedy(read, path, {(unitig.substr(unitig.size()-k+1,k-1)),overlap.second+(unitig.size()-k+1)},errors);
			//~ }else if(miss==minMiss){
				//~ equal=true;
			//~ }else if(miss<minMiss){
				//~ nextOverlap=((unitig.substr(unitig.size()-k+1,k-1)));
				//~ minMiss=miss;
				//~ indiceMinMiss=i;
				//~ nextUnitig=unitig;
				//~ equal=ended=false;
			//~ }
		//~ }
	//~ }
	//~ if(minMiss<=errors and not equal){
		//~ path.push_back(rangeUnitigs[indiceMinMiss].second);
		//~ if(ended){return minMiss;}
		//~ return minMiss+mapOnRightEndGreedy(read, path, {nextOverlap,overlap.second+(nextUnitig.size()-k+1)},errors-minMiss);
	//~ }
	//~ if(equal){
		//~ return 1000;
	//~ }
	//~ return minMiss;
//~ }



//~ vector<int> Aligner::inclued(vector<int>& v1, vector<int>& v2){
	//~ if(v1.size()<v2.size()){
		//~ for(uint i(1);i<v1.size();++i){
			//~ bool found (false);
			//~ for(uint j(1);j<v2.size();++j){
				//~ if(v1[i]==v2[j] or v1[i]==-v2[j]){
					//~ found=true;
					//~ break;
				//~ }
			//~ }
			//~ if(not found){
			//~ }
		//~ }
		//~ return v1;
	//~ }else{
		//~ for(uint i(1);i<v2.size();++i){
			//~ bool found (false);
			//~ for(uint j(1);j<v1.size();++j){
				//~ if(v2[i]==v1[j] or v2[i]==-v1[j]){
					//~ found=true;
					//~ break;
				//~ }
			//~ }
			//~ if(not found){
				//~ return {};
			//~ }
		//~ }
		//~ return v2;
	//~ }

	//~ return {};
//~ }



void get_consensus( string& consensus,const string& newread,const string& actualread, bool print=false){
	//~ cout<<consensus<<"\n"<<newread<<"\n"<<actualread<<"\n\n";
	for(uint i(0);i<consensus.size();++i){
		if(consensus[i]!=newread[i]){
			consensus[i]=(actualread[i]);
			if(print){
				cout<<'M';
			}
		}else{
			if(print){
				cout<<' ';
			}
		}
	}
}



//OK
void Aligner::alignReadOpti(const string& read, alignment& al, bool perfect=false){
	cout<<"************************************************************************************************"<<endl;
	al.path={};
	alignment best_al;
	int max_score((int)read.size()-5*errorsMax);
	//~ cout<<"ARO"<<max_score<<endl;
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	if(listAnchors.empty()){
		++noOverlapRead;
		++notAligned;
		return;
	}
	string superRead,superReadMem;
	random_shuffle ( listAnchors.begin(), listAnchors.end() );
	bool found(false);
	for(uint i(0);i<listAnchors.size();++i){
		al.path={};
		if(stringMode){
			//TODO STR
			//~ path=alignReadGreedyAnchorsstr(read,max_score,listAnchors[i],score);
		}else{
			al=alignReadGreedyAnchors(read,max_score,listAnchors[i]);
			if(al.score==read.size()){
				++alignedRead;
				return;
			}
		}
		//~ alignment_clean(al);
		if(al.score>=max_score){
			if(false and al.score==max_score and noMultiMapping){
				if(al.path!=best_al.path){
					found=true;
				}
			}else if (al.score>max_score){
				max_score=al.score;
				best_al=al;
				found=false;
			}
		}
	}
	if(found){
		al.path={};
		++notAligned;
		return;
	}
	al=best_al;
	if(not al.path.empty()){
		++alignedRead;
	}else{
		++notAligned;
	}
}





string Aligner::alignReadOpti_correction(const string& read, alignment& al){
	al.path={};
	alignment best_al;
	int max_score(read.size()-5*errorsMax-1);
	string consensus,new_correction;
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	if(listAnchors.empty()){
		++noOverlapRead;
		++notAligned;
		return read;
	}
	string superRead,superReadMem;
	random_shuffle ( listAnchors.begin(), listAnchors.end() );
	bool found(false);
	for(uint i(0);i<listAnchors.size();++i){
		al.path={};
		if(stringMode){
			//TODO STR
			//~ path=alignReadGreedyAnchorsstr(read,max_score,listAnchors[i],score);
		}else{
			al=alignReadGreedyAnchors(read,max_score,listAnchors[i]);
			if(al.score==read.size()){
				++alignedRead;
				return read;
			}
		}
		alignment_clean(al);
		if(al.score>=max_score){
			if(al.score==max_score and noMultiMapping){
				new_correction=get_corrected_read(al,read);

				if(consensus!=""){
					get_consensus(consensus,new_correction,read);
					if(consensus==read){
						return read;
					}else{
					}
				}else{
					consensus=new_correction;
				}
			}else if (al.score>max_score){
				max_score=al.score;
				best_al=al;
				consensus=get_corrected_read(al,read);
			}
		}
	}
	al=best_al;
	if(not al.path.empty()){
		++alignedRead;
		return consensus;
	}else{
		++notAligned;
		return read;
	}
}


string Aligner::alignReadOpti_correction2(const string& read, alignment& al){
	al.path={};
	alignment best_al;
	int max_score(read.size()-5*errorsMax-1);
	string consensus,new_correction;
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	if(listAnchors.empty()){
		++noOverlapRead;
		++notAligned;
		return read;
	}
	string superRead,superReadMem;
	random_shuffle ( listAnchors.begin(), listAnchors.end() );
	bool found(false);
	for(uint i(0);i<listAnchors.size();++i){
		al.path={};
		if(stringMode){
			//TODO STR
			//~ path=alignReadGreedyAnchorsstr(read,max_score,listAnchors[i],score);
		}else{
			al=alignReadGreedyAnchors(read,max_score,listAnchors[i]);
			if(al.score==read.size()){
				++alignedRead;
				return read;
			}
		}
		alignment_clean(al);
		if(al.score>=max_score){
			if(al.score==max_score and noMultiMapping){
				new_correction=get_corrected_read(al,read);


			}else if (al.score>max_score){
				max_score=al.score;
				best_al=al;
				consensus=get_corrected_read(al,read);
			}
		}
	}
	al=best_al;
	if(not al.path.empty()){
		++alignedRead;
		return consensus;
	}else{
		++notAligned;
		return read;
	}
}






//TODO URGENT
//~ void Aligner::alignReadOpti_correction(const string& read, string& corrected_seq){
	//~ exit(0);
	//~ vector<int> path={};
	//~ uint errors(0);
	//~ vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	//~ vector<int> errorsFromPreviousMapping(listAnchors.size(),-1);
	//~ if(listAnchors.empty()){
		//~ ++noOverlapRead;
		//~ return;
	//~ }
	//~ string superRead,superReadMem;
	//~ random_shuffle ( listAnchors.begin(), listAnchors.end() );
	//~ while(errors<=errorsMax){
		//~ bool found(false);
		//~ int errorInMapping(0);
		//~ for(uint i(0);i<listAnchors.size();++i){
			//~ path={};
			//~ errorInMapping=0;
			//~ if(errorsFromPreviousMapping[i]<=(int)errors){
				//~ if(stringMode){
					//path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i],errorInMapping);
				//~ }else{
					//~ path=alignReadGreedyAnchors(read,read.size(),listAnchors[i],errorInMapping);
				//~ }
				//~ errorsFromPreviousMapping[i]=errorInMapping;
				//~ //MAPPING IS FOUND
				//~ al=alignment_clean(al,read.size());
				//~ if(not path.empty()){
					//~ if(found){
						//~ superRead=(recoverSuperReadsCor(path,read.size()));
						//~ consensus(corrected_seq,superRead,read);
					//~ }else{
						//~ found=true;
						//~ superRead=(recoverSuperReadsCor(path,read.size()));
						//~ corrected_seq=superRead;
					//~ }
				//~ }
			//~ }
		//~ }
		//~ ++errors;
		//~ if(found){
			//~ ++alignedRead;
			//~ return;
		//~ }
	//~ }
	//~ if(not path.empty()){
		//~ ++alignedRead;
	//~ }
//~ }




//~ void Aligner::alignReadAllOpti(const string& read, vector<vector<int>>& pathVector){
	//~ pathVector={};
	//~ vector<int> path;
	//~ uint errors(0);
	//~ vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	//~ if(listAnchors.empty()){++noOverlapRead;return;}
	//~ while(errors<=errorsMax){
		//~ bool found(false);
		//~ int errorInMapping(0);
		//~ for(uint i(0);i<listAnchors.size();++i){
			//~ path={};
			//~ if(stringMode){
				//~ path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i],errorInMapping);
			//~ }else{
				//~ path=alignReadGreedyAnchors(read,errors,listAnchors[i],errorInMapping);
			//~ }
			//~ //MAPPING IS FOUND
			//~ if(not path.empty()){
				//~ pathVector.push_back(path);
				//~ found=true;
			//~ }
		//~ }
		//~ ++errors;
		//~ if(found){
			//~ ++alignedRead;
			//~ return;
		//~ }
	//~ }
//~ }



//~ void Aligner::alignReadAll(const string& read, vector<vector<int>>& pathVector){
	//~ pathVector={};
	//~ vector<int> path;
	//~ uint errors(0);
	//~ vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	//~ if(listAnchors.empty()){++noOverlapRead;return;}
	//~ while(errors<=errorsMax){
		//~ for(uint i(0);i<listAnchors.size();++i){
			//~ path={};
			//~ uint errorInMapping(0);
			//~ if(stringMode){
				//~ path=alignReadGreedyAnchorsstr(read,errors,listAnchors[i],errorInMapping);
			//~ }else{
				//~ path=alignReadGreedyAnchors(read,errors,listAnchors[i],errorInMapping);
			//~ }
			//~ //MAPPING IS FOUND
			//~ if(not path.empty()){
				//~ pathVector.push_back(path);
			//~ }
		//~ }
		//~ ++errors;
	//~ }
//~ }



void Aligner::alignReadFrom(const string& read, vector<int>& path, int unumber){
	vector<pair<pair<uint,uint>,uint>> listAnchors(getNAnchors(read,tryNumber));
	for(uint i(0);i<listAnchors.size();++i){
		auto  anchor=listAnchors[i];
		if(((anchor.first.first==unumber or anchor.first.first==-unumber) and anchor.second==0 and anchor.first.second==0)){
			path={};
			int errorInMapping(0);
			if(stringMode){
				//~ path=alignReadGreedyAnchorsstr(read,0,anchor,errorInMapping);
				//~ cout<<"NOPE"<<endl;exit(0);
			}else{
				alignment al=alignReadGreedyAnchors(read,read.size(),anchor);
				path=al.path;
			}
			if(not path.empty()){
				return;
			}
		}
	}
}



string Aligner::get_corrected_read(const alignment& al, const string& read){
	//~ cout<<al.score<<endl;
	if(al.path.empty()){
		cout<<"EMPTY"<<endl;
		cin.get();
		return read;
	}
	if(al.unitig_number_anchored>=al.path.size()){
		//~ cout<<"PB";
		cout<<"DAFUCK"<<endl;
		cin.get();
		return read;
	}
	//~ cout<<al.path.size()<<" "<<al.unitig_number_anchored<<endl;
	//~ cout<<1<<endl;

	//~ cout<<al.position_anchors_in_unitig<<" "<<al.position_anchors_in_read<<endl;
	string consensus(getUnitig(al.path[0])),unitig,inter;
	int position_start(0),nucleotides_in_path(0);
	if(al.unitig_number_anchored==0){
		//~ cout<<"found"<<endl;
		position_start=nucleotides_in_path+al.position_anchors_in_unitig-al.position_anchors_in_read;
	}
	nucleotides_in_path+=consensus.size()-k+1;
	//~ cout<<2<<endl;
	for(uint i(1); i<al.path.size(); ++i){
		//~ cout<<3<<endl;
		unitig=(getUnitig(al.path[i]));
		inter=(compactionEndNoRC(consensus, unitig, k-1));

		//~ cout<<"NIP"<<nucleotides_in_path<<endl;
		if(al.unitig_number_anchored==i){
			//~ cout<<"found"<<endl;
			//~ cout<<nucleotides_in_path<<" "<<al.position_anchors_in_read<<endl;;
			position_start=nucleotides_in_path+al.position_anchors_in_unitig-al.position_anchors_in_read;
			//~ cout<<position_start<<endl;
		}
		nucleotides_in_path+=unitig.size()-k+1;
		if(inter.empty()){
			cout<<i<<endl;cout<<"bug compaction super reads"<<endl;return {};
		}else{
			consensus=inter;
		}
	}
	//~ cout<<position_start<<endl;
	//~ cout<<consensus<<endl;
	if(position_start>=0){
		//~ cout<<"substr"<<endl;
		consensus=consensus.substr(position_start, read.size());
	}else{
		//~ cout<<"add prefix"<<endl;
		consensus=read.substr(0,-position_start)+consensus;
		consensus.substr(0,read.size());
	}
	if(consensus.size()<read.size()){
		//~ cout<<"add suffix"<<endl;
		//~ cout<<consensus<<endl;cin.get();
		consensus+=read.substr(consensus.size());
		//~ cout<<consensus<<endl;
	}

	//~ cin.get();
	//~ if((int)missmatchNumber(read,consensus,1000)<read.size()-25){
		//~ cout<<consensus<<endl;
		//~ cout<<missmatchNumber(read,consensus,1000)<<endl;
		//~ cout<<read<<endl;
		//~ cout<<consensus<<endl;
		//~ cin.get();
	//~ }
	int new_score(missmatchNumber(read,consensus,1000));
	//~ cout<<new_score<<" "<<al.score<<endl;
	if(al.score!=new_score){
		cout<<"WTFFFF"<<endl;
		cin.get();
	}
	return consensus;
}




void Aligner::alignPartGreedy(uint indiceThread){
	vector<pair<string,string>> multiread;
	//~ vector<uNumber> path,path2;
	alignment al,al2;
	vector<vector<int>> pathVector;
	string read,read2,header,header2,corrected,superRead,superRead2,suffix_stop,toWrite,align,consensus;
	vector<string> toWriteComp(nbBuckets);
	pair<string,string> superpath;
	while(not readFile->eof()){
	//~ while(not feof(readFileF)){
		toWrite="";
		if(keepOrder){
			while(threadToRead!=indiceThread){
				this_thread::sleep_for (chrono::microseconds(1));
			}
			readMutex.lock();
			getReads(multiread,100000);
			threadToRead=(threadToRead+1)%coreNumber;
			readMutex.unlock();
		}else{
			readMutex.lock();
			getReads(multiread,100000);
			readMutex.unlock();
		}
		if(pairedMode){
			//PAIRED READS
			for(uint i(0);i+1<multiread.size();i+=2){
				header=multiread[i].first;
				read=multiread[i].second;
				header2=multiread[i+1].first;
				read2=multiread[i+1].second;
				readNumber+=2;
				bool rc(false), noOverlap(false), overlapFound(false);
				alignReadOpti(read,al);
				alignReadOpti(read2,al2);
				if(al.path.empty()){
					++notAligned;
				}else{
					//~ path=cleanSR(path,read.size());
				}
				if(al2.path.empty()){
					++notAligned;
				}else{
					//~ path2=cleanSR(path2,read2.size());
				}
				auto superpath=(recoverSuperReadsPairedNoStr(al.path,al2.path));
				if(correctionMode){
					//~ cout<<1<<endl;
					if(al.path.empty()){
						if(al2.path.empty()){
							//~ cout<<2<<endl;
							toWrite+=header+'\n'+read+'\n';
							toWrite+=header2+'\n'+read2+'\n';
							continue;
						}else{
							//~ cout<<3<<endl;
							superRead2=(recoverSuperReadsCor(al2.path,read2.size()));
							toWrite+=header+'\n'+read+'\n';
							toWrite+=header2+'\n'+superRead2+'\n';
							continue;
						}
					}else{
						//~ cout<<4<<endl;
						if(al2.path.empty()){
							//~ cout<<5<<endl;
							superRead=(recoverSuperReadsCor(al.path,read.size()));
							toWrite+=header+'\n'+superRead+'\n';
							toWrite+=header2+'\n'+read2+'\n';
							continue;
						}
					}
					//BOTH ARE LAIGNED AT THIS POINT
					auto superpath_numbers=(recoverSuperReadsPaired_numbers(al.path,al2.path));
					if(not superpath_numbers.first.empty()){
						if(not superpath_numbers.second.empty()){
							//~ cout<<6<<endl;
							//NO MERGE
							superRead2=(recoverSuperReadsCor(al2.path,read2.size()));
							superRead=(recoverSuperReadsCor(al.path,read.size()));
							toWrite+=header+'\n'+superRead+'\n'+header2+'\n'+superRead2+'\n';
						}else{
							//~ cout<<7<<endl;
							//MERGE
							superRead=(recoverSuperReadsCor(al.path));
							toWrite+=header+'\n'+superRead+'\n';
						}
					}else{
						cout <<"should not happend"<<endl;
					}
					continue;
				}
				if(superpath.first!=""){
					if(superpath.second!=""){
						if(headerNeeded){
							toWrite+=header+'\n'+superpath.first+'\n'+header2+'\n'+superpath.second+'\n';
						}else{
							toWrite+=superpath.first+'\n'+superpath.second+'\n';
							//~ toWrite+=superpath.first+superpath.second;
						}
					}else{
						if(headerNeeded){
							toWrite+=header+'\n'+superpath.first+'\n';
						}else{
							toWrite+=superpath.first+'\n';
							//~ toWrite+=superpath.first;
						}
					}
				}else{
					if(superpath.second!=""){
						if(headerNeeded){
							toWrite+=header2+'\n'+superpath.second+'\n';
							//~ toWrite+=header2+superpath.second;
						}else{
							toWrite+=superpath.second+'\n';
							//~ toWrite+=superpath.second;
						}
					}
				}
			}
		}else{
			//UNPAIRED READS
			for(uint i(0);i<multiread.size();i++){
				header=multiread[i].first;
				read=multiread[i].second;
				++readNumber;
				if(correctionMode){
					//~ consensus="";
					consensus=alignReadOpti_correction(read,al);
					//~ cout<<"al"<<endl;
					//~ consensus=get_corrected_read(al,read);
					//~ cout<<"corr"<<endl;
					//~ alignReadOpti_correction(read,consensus);
					if(consensus.empty()){
						toWrite+=header+'\n'+read+'\n';
						//~ ++notAligned;continue;
					}else{
						toWrite+=header+'\n'+consensus+'\n';
						continue;
					}
				}else if(uniqueOptimalMappingMode){
					alignReadOpti(read,al);
					if(al.path.empty()){
						//~ cerr<<header+'\n'+read+'\n';
						if(correctionMode){
							toWrite+=header+'\n'+read+'\n';
						}
						++notAligned;
						continue;
					}
					//REGULAR MODE
					//~ path=cleanSR(path,read.size());
					//~ printPath(path);
					//~ path=getcleanPaths(path,false,true);
					//~ printPath(path);
					superRead=(recoverSuperReadsNoStr(al.path,0));
					//~ cout<<superRead<<endl;cin.get();
					if(superRead!=""){
						if(printAlignment){
							toWrite+=header+'\n'+superRead+'\n';
							toWrite+=read+'\n'+recoverSuperReadsCor(al.path,read.size())+'\n';
						}else{
							if(headerNeeded){
								toWrite+=header+'\n'+superRead+'\n';
							}else{
								toWrite+=superRead+'\n';
							}
						}
					}
				}
			}
		}
		if(keepOrder){
			while(threadToPrint!=indiceThread){this_thread::sleep_for (chrono::microseconds(1));}
			pathMutex.lock();
			{
				//~ fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
				if(compression){
					pathCompressed->write(toWrite.c_str(),toWrite.size());
				}else{
					fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
				}
				threadToPrint=(threadToPrint+1)%coreNumber;
			}
			pathMutex.unlock();
		}else{
			if(not compressionMode){
				pathMutex.lock();
				{
					if(compression){
						pathCompressed->write(toWrite.c_str(),toWrite.size());
					}else{
						fwrite((toWrite).c_str(), sizeof(char), toWrite.size(), pathFilef);
					}
				}
				pathMutex.unlock();
			}else{
				for(uint i(0);i<nbBuckets;++i){
					pathMutexComp[i].lock();
					{
						fwrite((toWriteComp[i]).c_str(), sizeof(char), toWriteComp[i].size(), pathFileComp[i]);
						toWriteComp[i]="";
					}
					pathMutexComp[i].unlock();
				}
			}
		}
		progressMutex.lock();
		if(++iterLoop%100==0){
			cout<<"Reads: "<<intToString(readNumber)<<endl;
			cout<<"Not anchored : "<<intToString(noOverlapRead)<<" Percent: "<<(100*float(noOverlapRead))/readNumber<<endl;
			cout<<"Anchored and aligned : "<<intToString(alignedRead)<<" Percent: "<<(100*float(alignedRead))/(alignedRead+notAligned)<<endl;
			cout<<"Not aligned: "<<intToString(notAligned)<<" Percent: "<<(100*float(notAligned))/(alignedRead+notAligned)<<endl;
			if(pairedMode){
				cout<<"Super reads: "<<intToString(superReads)<<endl;
				cout<<"Failed super reads compaction: "<<intToString(notCompatedSR)<<endl;
			}
			cout<<endl;
		}else{
		}
		progressMutex.unlock();
	}
}



