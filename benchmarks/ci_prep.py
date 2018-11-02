import os
import shutil
import platform
import numpy as np
import pandas as pd
import platform
import pyemu

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"linux")
elif "darwin" in platform.platform().lower():
    bin_path = os.path.join(bin_path,"mac")
else:
    bin_path = os.path.join(bin_path,"win")

bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path


bin_path = os.path.join("..","..","..","bin")
    
if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "win", "pestpp-ies.exe")
elif "darwin" in platform.platform().lower():
    exe_path = os.path.join(bin_path,  "mac", "pestpp-ies")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-ies")

noptmax = 4
num_reals = 20
port = 4020

def setup_suite_dir(model_d):
    pyemu.Ensemble.reseed()
    base_d = os.path.join(model_d, "template")
    new_d = os.path.join(model_d, "test_template")
    if os.path.exists(new_d):
        shutil.rmtree(new_d)
    shutil.copytree(base_d, new_d)
    print(platform.platform().lower())
    local=True
    if "linux" in platform.platform().lower() and "10par" in model_d:
        #print("travis_prep")
        #prep_for_travis(model_d)
        local=False
    pst = pyemu.Pst(os.path.join(new_d, "pest.pst"))
    print(pst.model_command)
    
    # set first par as fixed
    #pst.parameter_data.loc[pst.par_names[0], "partrans"] = "fixed"

    pst.observation_data.loc[pst.nnz_obs_names,"weight"] = 1.0

    # set noptmax
    pst.control_data.noptmax = noptmax

    # wipe all pestpp options
    pst.pestpp_options = {}
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["lambda_scale_fac"] = 1.0
    pst.pestpp_options["ies_lambda_mults"] = 1.0
    # write a generic 2D cov
    if os.path.exists(os.path.join(new_d,"prior.jcb")):
        cov = pyemu.Cov.from_binary(os.path.join(new_d,"prior.jcb"))
        #cov.to_ascii(os.path.join(new_d,"prior.cov"))
    elif os.path.exists(os.path.join(new_d, "prior.cov")):
        cov = pyemu.Cov.from_ascii(os.path.join(new_d, "prior.cov"))
    else:
        cov = pyemu.Cov.from_parameter_data(pst)
        cov = pyemu.Cov(cov.as_2d, names=cov.row_names)
        #cov.to_ascii(os.path.join(new_d, "prior.cov"))
        cov.to_binary(os.path.join(new_d, "prior.jcb"))

    # draw some ensembles
    idx = [i for i in range(num_reals)]
    idx[-1] = "base"
    pe = pyemu.ParameterEnsemble.from_gaussian_draw(pst, cov=cov, num_reals=num_reals,
                                                    use_homegrown=True,group_chunks=True)
    pe.index = idx
    pe.to_csv(os.path.join(new_d, "par.csv"))
    pe.to_binary(os.path.join(new_d, "par.jcb"))
    pe.to_csv(os.path.join(new_d, "sweep_in.csv"))
    pe.loc[:, pst.adj_par_names].to_csv(os.path.join(new_d, "par_some.csv"))
    pe.iloc[:-3, :].to_csv(os.path.join(new_d, "restart_failed_par.csv"))
    oe = pyemu.ObservationEnsemble.from_id_gaussian_draw(pst, num_reals=num_reals)
    oe.index = idx
    oe.to_csv(os.path.join(new_d, "obs.csv"))
    oe.iloc[:-3, :].to_csv(os.path.join(new_d, "restart_failed_base_obs.csv"))
    oe.to_binary(os.path.join(new_d, "obs.jcb"))

    pst.write(os.path.join(new_d, "pest.pst"))
    # run sweep
    m_d = os.path.join(model_d,"master_sweep1")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_slaves(new_d, exe_path.replace("-ies","-swp"), "pest.pst", 5, master_dir=m_d,
                           slave_root=model_d,local=local,port=port)
    df = pd.read_csv(os.path.join(m_d, "sweep_out.csv"),index_col=0)
    
    m_d = os.path.join(model_d,"master_pestpp")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_slaves(new_d, exe_path.replace("-ies",""), "pest.pst", 10, master_dir=m_d,
                           slave_root=model_d,local=local,port=port)

    m_d = os.path.join(model_d,"master_pestpp-ies")
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    pyemu.os_utils.start_slaves(new_d, exe_path, "pest.pst", 10, master_dir=m_d,
                           slave_root=model_d,local=local,port=port)



    

if __name__ == "__main__":
	setup_suite_dir("ies_10par_xsec")
