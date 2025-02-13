/**
 * @file
 * @brief Others
 */
 
#include "HEP.h"
#include <cmath>

namespace HepLib {
    
    namespace {
        class Family {
            public:
                int pn;
                lst intg;
                Family(int pn_) : pn(pn_) { }
            };
    }
    
    void Export2AMFlow(const ex & expr, const string & dir, const lst & loop, const lst & leg, const lst & replacement, const lst & numeric) {
        static string wlo = R"EOF(
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;
<<AMFlow`AMFlow`
SetReductionOptions["IBPReducer"->"FiniteFlow+LiteRed"];
ID = ToExpression[$ScriptCommandLine[[2]]];
PropsInts=Get[FileNameJoin[{current,ToString[ID]<>".m"}]];

AMFlowInfo["Family"] = Symbol["I"<>ToString[ID]];
AMFlowInfo["Loop"] = <<Loop>>;
AMFlowInfo["Leg"] = <<Leg>>;
AMFlowInfo["Conservation"] = { };
AMFlowInfo["Replacement"] = <<Replacement>>;
AMFlowInfo["Propagator"] = PropsInts[[1]];
AMFlowInfo["Numeric"] = <<Numeric>>;
AMFlowInfo["NThread"] = 6;

integrals = Map[(j[Symbol["I"<>ToString[ID]], Sequence@@#])&, PropsInts[[2]]];
precision = 30;
epsorder = 2*Length@AMFlowInfo["Loop"]-PropsInts[[3]];
sol = SolveIntegrals[integrals, precision, epsorder];
Put[sol, FileNameJoin[{current, "sol"<>ToString[ID]<>".m"}]];

Quit[];
)EOF";

        string sLoop = ex2str(loop);
        string_replace_all(sLoop, "==", "->");
        string sLeg = ex2str(leg);
        string_replace_all(sLeg, "==", "->");
        string sReplacement = ex2str(replacement);
        string_replace_all(sReplacement, "==", "->");
        string sNumeric = ex2str(numeric);
        string_replace_all(sNumeric, "==", "->");
        string wl = wlo;
        string_replace_all(wl, "<<Loop>>", sLoop);
        string_replace_all(wl, "<<Leg>>", sLeg);
        string_replace_all(wl, "<<Replacement>>", sReplacement);
        string_replace_all(wl, "<<Numeric>>", sNumeric);

        system(("mkdir -p "+dir).c_str());
        str2file(wl, dir+"/run.wl");
        
        auto res = collect_ex(expr, F(w1,w2));
        
        exset fs;
        find(expr, F(w1,w2), fs);
        map<ex,Family,ex_is_less> p2f;
        exmap f2f;
        int pn = 0;
        for(auto f : fs) {
            auto itr = p2f.find(f.op(0));
            if(itr == p2f.end()) {
                itr = p2f.insert({f.op(0), Family(pn)}).first;
                pn++;
            }
            itr->second.intg.append(f.op(1));
            f2f[f] = F(itr->second.pn, f.op(1));
        }
        res = res.subs(f2f);
        
        ex2file(res, dir+"/res.txt");
        
        string sres = R"EOF(
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;
<<HepLib`
Module[{sols, res},
sols = Flatten@Table[Import[FileNameJoin[{current,"sol"<>ToString[i]<>".m"}]], {i,0,<<NN>>}];
res = C2M@Import[FileNameJoin[{current, "res.txt"}]];
res = res/.{F[p_, n_] :> j[Symbol["I"<>ToString[p]], Sequence@@n]} /. sols;
res//Print
]
)EOF";
        
        string_replace_all(sres, "<<NN>>", to_string(pn-1));
        str2file(sres, dir+"/res.m");
        
        for(auto kv : p2f) {
            lst prop = ex_to<lst>(kv.first);
            auto & f = kv.second;
            for(int i=0; i<prop.nops(); i++) {
                bool ok = false;
                for(int j=1; j<=4; j++) {
                    Symbol qj("q"+to_string(j));
                    if(prop.op(i).has(qj*qj)) {
                        ok = true;
                        break;
                    }
                }
                if(!ok) {
                    cout << "Propagators: " << prop.op(i) << " @ pos: " << i << endl;
                    cout << "Integrals: " << f.intg << endl;
                    for(auto item : f.intg) {
                        if(!item.op(i).is_zero()) {
                            ok = true;
                            break;
                        }
                    }
                    if(ok) {
                        cout << "Error: index is NOT 0, abort!" << endl;
                        abort();
                    } else {
                        cout << WHITE << "ok: " << RESET << prop.op(i);
                        prop.let_op(i) = prop.op(i).subs(w1*w2 == pow(w1+11*w2,2));
                        cout << " -> " << prop.op(i) << endl;
                        cout << endl;
                    }
                }
            }
            
            int order = 0;
            for(auto ns : f.intg) {
                ex cc = res.coeff(F(f.pn, ns));
                cc = cc.subs(d==4-2*ep);
                cc = series_ex(cc, ep, 0);
                if(cc.is_zero()) cc = 1;
                auto ldeg = cc.ldegree(ep);
                if(order>ldeg) order = ldeg;
            }
            
            ex2file(lst{ prop, f.intg, order }, dir+"/"+to_string(f.pn)+".m");
        }
                
    }

}
