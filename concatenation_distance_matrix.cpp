#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>

using namespace std;

int dp[20000][20000];
int mat[7][5][5];
const int inf=1e9;
string con[5];

map<pair<string,string>, pair<int,int> > m;  

int align(string s, string t){
	int match=0,mismatch=1,gap=1,ret=inf;
	for(int i=0;i<=s.size();i++)for(int j=0;j<=t.size();j++)dp[i][j]=inf;
	for(int i=0;i<=s.size();i++)dp[i][0]=i;
	for(int j=1;j<=t.size();j++)dp[0][j]=j;
	for(int i=1;i<=s.size();i++){
		for(int j=1;j<=t.size();j++){
			dp[i][j]=min(dp[i-1][j],dp[i][j-1])+gap;
			if(s[i-1]==t[j-1])dp[i][j]=min(dp[i][j],dp[i-1][j-1]+match);
			else dp[i][j]=min(dp[i][j],dp[i-1][j-1]+mismatch);
		}
	}
	return dp[s.size()][t.size()];
}

string names[]={
	"Zaire", "TaiForest", "Sudan", "Reston", "Bundibugyo"
};
string Gene[]={
	"NP", "VP35", "VP40", "GP", "VP30", "VP24", "L"
};

int main(){
	ifstream cin("output/gene_alignment");
	string x,y;
	int lo,hi;
	while(cin>>x>>y>>lo>>hi){
		m[make_pair(x,y)]=make_pair(lo,hi);
	}
	for(int i=0;i<5;i++){
		string name="resources/"+names[i]+"_genome.fasta";
		ifstream cin(name);
		string t="";
		string tmp;
		cin>>tmp;
		string last=tmp.substr(1,tmp.size()-1);
		getline(cin,tmp);
		while(getline(cin,tmp)){
			tmp=tmp.substr(0,tmp.size()-1);
			t+=tmp;
		}
		string left=t;
		for(int k=0;k<7;k++){
			pair<int,int> range_left=m[make_pair(names[i],Gene[k])];
			string new_left=left.substr(range_left.first,range_left.second-range_left.first);
			con[i]+=new_left;
		}
	}
	ofstream cout("output/concatenation.csv");
	for(int i=0;i<5;i++){
		cout<<names[i];
		if(i!=4)cout<<",";
	}
	cout<<endl;
	for(int i=0;i<5;i++){
		for(int j=0;j<5;j++){
			cout<<align(con[i],con[j]);
			if(j!=4)cout<<",";
		}
		cout<<endl;
	}
	return 0;
}
