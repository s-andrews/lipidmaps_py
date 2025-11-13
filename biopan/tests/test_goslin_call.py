from biopan.src.models.goslin import Goslin
class Sample:
    def __init__(self, name):
        self.sample_name = name

samps = [Sample("Cholesterol"), Sample("PC(16:0/18:1)")]
res = Goslin.annotate_samples(samps)
print('lookup:', res)
for s in samps:
    print(s.sample_name, '->', getattr(s, 'goslin_result', None))