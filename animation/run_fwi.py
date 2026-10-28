"""Run the 3 -> 5 -> 8 Hz Marmousi FWI and save the model after every L-BFGS update
(plus an 8 Hz RTM image) to fwi_frames.npz for the time-lapse video."""
import sys, os, time, numpy as np
DEVITO_REPO = os.environ.get('DEVITO_REPO', './devito')         # cloned devitocodes/devito (v4.8.23)
DATA_PATH = os.environ.get('DATA_PATH', './devito-data')        # cloned devitocodes/data
sys.path.insert(0, DEVITO_REPO)
from devito import configuration, Function, norm
configuration['log-level'] = 'ERROR'
from examples.seismic import SeismicModel, AcquisitionGeometry, Receiver
from examples.seismic.acoustic import AcousticWaveSolver
from scipy.ndimage import gaussian_filter, zoom, laplace
from scipy.optimize import minimize
vm = np.fromfile(f'{DATA_PATH}/Simple2D/vp_marmousi_bi', dtype='float32').reshape(1601, 401)
vt = zoom(vm[600:1200, :240], 7.5/25, order=1).astype(np.float32)
dx, nbl, so, DT = 25., 40, 8, 2.3
water = vt <= 1.51; upd = ~water
nx, nz = vt.shape; Lx = (nx-1)*dx
def make_model(v, grid=None):
    return SeismicModel(vp=v, origin=(0.,0.), spacing=(dx,dx), shape=v.shape, space_order=so, nbl=nbl, bcs='damp', grid=grid, dt=DT)
model = make_model(vt)
model0 = make_model(np.where(water, 1.5, gaussian_filter(vt, 6)).astype(np.float32), grid=model.grid)
v_init = np.asarray(model0.vp.data)[nbl:-nbl, nbl:-nbl].copy()
nrec, nshots = 181, 16
rec_xy = np.zeros((nrec,2)); rec_xy[:,0] = np.linspace(0, Lx, nrec); rec_xy[:,1] = 12.5
src_x = np.linspace(100, Lx-100, nshots)
def geometry(f): return AcquisitionGeometry(model, rec_xy, np.array([[Lx/2, 12.5]]), 0., 3000., f0=f/1000, src_type='Ricker')
def rel_err(v): return np.linalg.norm((v-vt)[upd])/np.linalg.norm(vt[upd])*100
class FWI:
    def __init__(s, f):
        s.geo = geometry(f); s.solver = AcousticWaveSolver(model, s.geo, space_order=so)
        s.res = Receiver(name='res', grid=model.grid, time_range=s.geo.time_axis, coordinates=s.geo.rec_positions)
        s.dobs = []
        for sx in src_x:
            s.geo.src_positions[0,0] = sx; s.dobs.append(np.array(s.solver.forward(vp=model.vp)[0].data))
    def __call__(s, v, illum=False):
        model0.update('vp', v); grad = Function(name='grad', grid=model.grid); J = 0.; E = 0.
        for sx, dob in zip(src_x, s.dobs):
            s.geo.src_positions[0,0] = sx
            dsyn, u0, _ = s.solver.forward(vp=model0.vp, save=True)
            s.res.data[:] = np.asarray(dsyn.data) - dob; J += .5*norm(s.res)**2
            s.solver.jacobian_adjoint(s.res, u0, vp=model0.vp, grad=grad)
            if illum: E = E + np.sum(np.asarray(u0.data)[:, nbl:-nbl, nbl:-nbl]**2, axis=0)
        gv = np.asarray(grad.data)[nbl:-nbl, nbl:-nbl].copy()*(-2.0/v**3)
        return (J, gv, E) if illum else (J, gv)
snaps = [v_init.copy()]; labels = [(0, 0)]; errs = [rel_err(v_init)]
v = v_init.copy()
for f0 in [3, 5, 8]:
    t0 = time.time(); fwi = FWI(f0); _, _, E = fwi(v, illum=True)
    P = 1.0/np.sqrt(E/E[upd].max()+0.05); P = P/P[upd].min(); p = P[upd]; vf = v.copy()
    def f(x):
        vf[upd] = x*p; J, gv = fwi(vf.astype(np.float32)); return J, (gv[upd]*p).astype(np.float64)
    k = [0]
    def cb(x):
        k[0] += 1; vv = v.copy(); vv[upd] = x*p; snaps.append(vv.astype(np.float32)); labels.append((f0, k[0])); errs.append(rel_err(vv))
        print(f0, k[0], f'{errs[-1]:.2f}%', flush=True)
    r = minimize(f, v[upd]/p, jac=True, method='L-BFGS-B', bounds=list(zip(1.5/p, 4.7/p)), callback=cb, options=dict(maxiter=12, maxfun=20))
    v = v.copy(); v[upd] = r.x*p; v = v.astype(np.float32)
    print(f'band {f0} done {time.time()-t0:.0f}s err {rel_err(v):.2f}%', flush=True)
# RTM image with final model, 8 Hz, remove direct wave by subtracting data from water-only model
geo = geometry(8); sol = AcousticWaveSolver(model, geo, space_order=so)
mw = make_model(np.full_like(vt, 1.5), grid=model.grid)
model0.update('vp', v); img = Function(name='img', grid=model.grid)
res = Receiver(name='r2', grid=model.grid, time_range=geo.time_axis, coordinates=geo.rec_positions)
for sx in src_x:
    geo.src_positions[0,0] = sx
    dob = np.array(sol.forward(vp=model.vp)[0].data); dw = np.array(sol.forward(vp=mw.vp)[0].data)
    _, u0, _ = sol.forward(vp=model0.vp, save=True)
    res.data[:] = dob - dw
    sol.jacobian_adjoint(res, u0, vp=model0.vp, grad=img)
im = -laplace(np.asarray(img.data)[nbl:-nbl, nbl:-nbl])
np.savez('fwi_frames.npz', snaps=np.array(snaps), labels=np.array(labels), errs=np.array(errs), vt=vt, v0=v_init, img=im, water=water)
print('saved', len(snaps))
