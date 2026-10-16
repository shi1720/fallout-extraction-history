"""Independent finite feasibility checks and source-arithmetic tests."""
import sys, unittest, math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.optimize import linprog
from analyze import read_soils, fit_models, proportional_interval, additive_interval, minimum_tolerance

def lp_feasible(x,y,e,kind):
    # Latent quantities are explicit variables. The oracle never calls interval formulas.
    n=len(x); A=[];b=[]
    if kind=='proportional':
        bounds=[(0,1)]+[(yy*(1-e),yy*(1+e)) for yy in y]
        for j,xx in enumerate(x,1):
            row=np.zeros(n+1);row[0]=xx*(1-e);row[j]=-1;A.append(row);b.append(0)
            row=np.zeros(n+1);row[0]=-xx*(1+e);row[j]=1;A.append(row);b.append(0)
    else:
        bounds=[(0,None)]+[(xx*(1-e),xx*(1+e)) for xx in x]
        for j,yy in enumerate(y,1):
            row=np.zeros(n+1);row[j]=1;row[0]=-1;A.append(row);b.append(yy*(1+e))
            row=-row;A.append(row);b.append(-yy*(1-e))
    result=linprog(np.zeros(n+1),A_ub=A,b_ub=b,bounds=bounds,method='highs')
    if result.status not in (0,2):raise AssertionError(result.message)
    return result.status==0

class ReanalysisTests(unittest.TestCase):
    def setUp(self):
        self.rows=read_soils();self.x=np.array([float(r['fusion_mCi_mi2']) for r in self.rows]);self.y=np.array([float(r['hcl_mCi_mi2']) for r in self.rows])
    def test_source_enumeration(self):
        self.assertEqual([int(r['site_id']) for r in self.rows],list(range(1,17)))
        self.assertEqual(sum(self.x),898);self.assertEqual(sum(self.y),797)
    def test_printed_ratios(self):
        for row,x,y in zip(self.rows,self.x,self.y):self.assertEqual(round(y/x*100),int(row['printed_percent']))
    def test_range_discrepancy_survives_rounding(self):
        self.assertEqual(int(self.rows[9]['printed_percent']),67)
        self.assertLess((14+.5)/(21-.5)*100,78)
    def test_regression_normal_equations(self):
        f=fit_models(self.x,self.y)['affine'];r=self.y-f['intercept']-f['slope']*self.x
        self.assertAlmostEqual(sum(r),0,places=9);self.assertAlmostEqual(sum(r*self.x),0,places=8)
        self.assertAlmostEqual(f['x_intercept'],2.502800839816,places=9)
    def test_rounding_can_explain_historical_sse(self):
        f=fit_models(self.x,self.y)['affine'];p=f['intercept']+f['slope']*self.x
        witness=p+math.sqrt(359/f['sse'])*(self.y-p)
        self.assertTrue(np.all(np.abs(witness-self.y)<.5))
        self.assertAlmostEqual(fit_models(self.x,witness)['affine']['sse'],359,places=8)
    def test_known_tolerance_boundaries(self):
        for m,fun in [('proportional',proportional_interval),('additive',additive_interval)]:
            t=minimum_tolerance(self.x,self.y,m)
            self.assertFalse(lp_feasible(self.x,self.y,t-1e-6,m));self.assertTrue(lp_feasible(self.x,self.y,t+1e-6,m))
            self.assertGreater(fun(self.x,self.y,t-1e-6)[0],fun(self.x,self.y,t-1e-6)[1])
    def test_ten_percent_compatibility(self):
        self.assertFalse(lp_feasible(self.x,self.y,.1,'proportional'));self.assertTrue(lp_feasible(self.x,self.y,.1,'additive'))
        self.assertTrue(lp_feasible(np.delete(self.x,9),np.delete(self.y,9),.1,'proportional'))
    def test_randomized_against_independent_lp(self):
        rng=np.random.default_rng(20261016)
        for case in range(300):
            n=int(rng.integers(2,12));x=rng.uniform(.5,100,n);y=rng.uniform(.5,100,n);e=float(rng.uniform(0,.8))
            for m,fun in [('proportional',proportional_interval),('additive',additive_interval)]:
                lo,hi=fun(x,y,e);self.assertEqual(lo<=hi,lp_feasible(x,y,e,m),(case,m))
                t=minimum_tolerance(x,y,m)
                self.assertFalse(lp_feasible(x,y,max(0,t-1e-5),m))
                self.assertTrue(lp_feasible(x,y,t+1e-5,m))
    def test_invalid_tolerance(self):
        for e in [-.1,1,2]:
            for f in [proportional_interval,additive_interval]:
                with self.assertRaises(ValueError):f(self.x,self.y,e)
if __name__=='__main__':unittest.main(verbosity=2)
