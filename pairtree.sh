module load pairtree/1.0.1

pat=$1

# run pairtree
pairtree --params ${pat}/mut_dyn/pairtree/json_${pat}_nc.json ${pat}/mut_dyn/pairtree/ssm_${pat}_nc.tsv ${pat}/mut_dyn/pairtree/${pat}_nc.npz
# visualize candidate trees and posterior
summposterior ${pat}/mut_dyn/pairtree/ssm_${pat}_nc.tsv ${pat}/mut_dyn/pairtree/json_${pat}_nc.json ${pat}/mut_dyn/pairtree/${pat}_nc.npz ${pat}/mut_dyn/pairtree/${pat}_nc_summposterior.html
# visualize highest posterior tree
plottree --tree-json ${pat}/mut_dyn/pairtree/plottree_${pat}_nc.json ${pat}/mut_dyn/pairtree/ssm_${pat}_nc.tsv ${pat}/mut_dyn/pairtree/json_${pat}_nc.json ${pat}/mut_dyn/pairtree/${pat}_nc.npz ${pat}/mut_dyn/pairtree/${pat}__nc_plottree.html
