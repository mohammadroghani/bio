#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

using namespace std;

int dp[20000][8000];
const int inf=1e9;


pair<int,int> align(string s, string t){
	int match=1,mismatch=-2,gap=-3,ret=-1;
	for(int i=0;i<=s.size();i++)for(int j=0;j<=t.size();j++)dp[i][j]=-inf;
	for(int i=0;i<=s.size();i++)dp[i][0]=0;
	for(int j=1;j<=t.size();j++)dp[0][j]=dp[0][j-1]+gap;
	for(int i=1;i<=s.size();i++){
		for(int j=1;j<=t.size();j++){
			dp[i][j]=max(dp[i-1][j],dp[i][j-1])+gap;
			if(s[i-1]==t[j-1])dp[i][j]=max(dp[i][j],dp[i-1][j-1]+match);
			else dp[i][j]=max(dp[i][j],dp[i-1][j-1]+mismatch);
		}
	}
	int mx=-inf;
	int x=-1,y=t.size();
	for(int i=0;i<=s.size();i++){
		if(dp[i][t.size()]>mx){
			mx=dp[i][t.size()];
			x=i;
		}
	}
	ret=x;
	while(y>0){
		if(x>0 && dp[x-1][y]+gap==dp[x][y]){
			x--;
		}
		else if(y>0 && dp[x][y-1]+gap==dp[x][y]){
			y--;
		}
		else if(x>0 && y>0 && s[x-1]==t[y-1] && dp[x][y]==dp[x-1][y-1]+match){
			x--;
			y--;
		}
		else if(x>0 && y>0 && s[x-1]!=t[y-1] && dp[x][y]==dp[x-1][y-1]+mismatch){
			x--;
			y--;
		}
	}
	return make_pair(x,ret);
}

string names[]={
	"Zaire", "TaiForest", "Sudan", "Reston", "Bundibugyo"
};

int main(){
	vector<pair<string,string> > v;
	string last="";
	string name="resources/Marburg_Genes.fasta";
	ifstream cin(name);
	ofstream cout("output/gene_alignment");
	string tmp,t;
	cin>>tmp;
	last=tmp.substr(1,tmp.size()-1);
	getline(cin,tmp);
	while(getline(cin,tmp)){
		if(tmp[0]=='>'){
			v.push_back(make_pair(last,t));
			t="";
			last=tmp.substr(1,tmp.size()-2);
			continue;
		}
		tmp=tmp.substr(0,tmp.size()-1);
		t+=tmp;
	}
	v.push_back(make_pair(last,t));
	for(int i=0;i<5;i++){
		for(int j=0;j<7;j++){
			name="resources/"+names[i]+"_genome.fasta";
			ifstream cin(name);
			t="";
			cin>>tmp;
			last=tmp.substr(1,tmp.size()-1);
			getline(cin,tmp);
			while(getline(cin,tmp)){
				tmp=tmp.substr(0,tmp.size()-1);
				t+=tmp;
			}
			pair<int,int> pos=align(t,v[j].second);
			cout<<names[i]<<" "<<v[j].first<<" "<<pos.first<<" "<<pos.second<<endl;
		}
	}
	return 0;
}
