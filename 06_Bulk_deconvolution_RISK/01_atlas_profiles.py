# Atlas UMI counts summed per (annotations_fine, donor)

import os
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scipy.io as sio

h5, out = os.environ["ATLAS_H5"], os.environ["OUT_DIR"]
os.makedirs(out, exist_ok=True)

a = ad.read_h5ad(h5)
obs = a.obs[["source", "participant_id", "patient", "Patient",
             "annotations_fine", "annotations_main"]].astype(str)
src = obs["source"]
donor = np.select([src == "PREDICT 5p", src == "PREDICT 3p", src == "Thomas"],
                  [obs["participant_id"], obs["patient"], obs["Patient"]], default="NA")
assert not np.isin(donor, ["NA", "nan"]).any(), "cells without donor id"
obs["donor"] = src + "|" + donor
print("cells:", len(obs), "| donors:", obs.donor.nunique(), flush=True)

X = a.raw.X.tocsr()
assert np.allclose(X.data[:100000], np.rint(X.data[:100000])), "raw.X is not integer counts"
keep = np.asarray(X.sum(0)).ravel() > 0
pd.Series(np.asarray(a.raw.var_names)[keep]).to_csv(f"{out}/genes.tsv", index=False, header=False)

key = obs["annotations_fine"] + "---" + obs["donor"]
codes, profiles = pd.factorize(key, sort=True)
G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                  shape=(len(profiles), len(codes)))
counts = (G @ X).tocsc()[:, keep]
counts.data = np.rint(counts.data)
sio.mmwrite(f"{out}/profiles_counts.mtx", counts.T.astype(np.int64), field="integer")

(obs.assign(profile=key)
    .groupby("profile")
    .agg(annotations_fine=("annotations_fine", "first"), main=("annotations_main", "first"),
         donor=("donor", "first"), source=("source", "first"), n_cells=("donor", "size"))
    .reindex(profiles)
    .to_csv(f"{out}/profiles.csv", index_label="profile"))
print("genes:", keep.sum(), "| profiles:", len(profiles))
