# getSignatures

```shell
ln -s ../SomatiCalling/results/ vcfs_raw
mkdir vcfs_PASS

for i in vcfs_raw/*norm.vcf.gz;
do
name=`basename $i|cut -d '.' -f 1`
if [ ! -f vcfs_PASS/${name}.PASS.vcf ];
then
    zgrep "PASS" $i | grep -v "#" > vcfs_PASS/${name}.PASS.vcf
fi
done
```

## SBS

```python
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import glob


Analyze.cosmic_fit(samples="vcfs_PASS",
                output = "SigProfilerAssignment_out/SBS,
                input_type = "vcf",
                context_type = "96",
                genome_build = "GRCh38",
                cosmic_version = 3.3,
                exclude_signature_subgroups = None,
                export_probabilities = True,
                export_probabilities_per_mutation = True,
                make_plots = True,
                sample_reconstruction_plots = True,
                verbose = True)
```