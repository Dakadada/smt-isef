"""Frozen, biologically typed nonstructural feature definitions."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class SpliceSite:
    position: int
    site_type: str
    endpoint_direction: int
    def __post_init__(self):
        if self.site_type not in ("donor", "acceptor") or self.endpoint_direction not in (-1, 1):
            raise ValueError("splice site needs donor/acceptor type and endpoint direction +/-1")


class InterventionInput:
    __slots__ = ("family", "study", "sample", "intervention", "event", "locus", "sequence", "footprint_start", "footprint_end", "splice_sites", "oligo_sequence", "chemistry", "dose", "footprint_kernel", "monovalent_salt_molar")
    def __init__(self, family, study, sample, intervention, event, locus, sequence, footprint_start, footprint_end, splice_sites, oligo_sequence, chemistry, dose=0.0, footprint_kernel=(), monovalent_salt_molar=0.15):
        self.family, self.study, self.sample, self.intervention = family, study, sample, intervention
        self.event, self.locus, self.sequence = event, locus, sequence
        self.footprint_start, self.footprint_end = int(footprint_start), int(footprint_end)
        self.splice_sites = tuple(splice_sites)
        if not self.splice_sites or not all(isinstance(x, SpliceSite) for x in self.splice_sites): raise ValueError("typed SpliceSite annotations are required")
        self.oligo_sequence, self.chemistry = oligo_sequence, chemistry
        self.dose, self.footprint_kernel, self.monovalent_salt_molar = float(dose), tuple(footprint_kernel), float(monovalent_salt_molar)


def _channels(backend, sequence):
    values = backend.score_channels(sequence)
    if not isinstance(values, dict) or set(values) != {"donor", "acceptor"}: raise ValueError("cis backend lacks typed donor/acceptor channels")
    if any(len(v) != len(sequence) for v in values.values()): raise ValueError("cis backend channel length mismatch")
    return values


def _event_score(channels, splice_sites):
    mapped = []
    for site in splice_sites:
        if not 0 <= site.position < len(channels[site.site_type]): raise ValueError("splice-site index outside sequence")
        mapped.append(site.endpoint_direction * channels[site.site_type][site.position])
    return sum(mapped), tuple(mapped)


def occlusion_importance(sequence, start, end, splice_sites, backend, kernel=()):
    """Frozen substitutions mapped channel-by-channel onto one endpoint direction."""
    seq = sequence.upper().replace("U", "T")
    if not (0 <= start < end <= len(seq)): raise ValueError("invalid half-open footprint")
    base_total, base_sites = _event_score(_channels(backend, seq), splice_sites)
    site_position_values = [[] for _ in splice_sites]
    for pos in range(start, end):
        per_site_changes = [[] for _ in splice_sites]
        for alt in "ACGT":
            if alt == seq[pos]: continue
            _, changed_sites = _event_score(_channels(backend, seq[:pos] + alt + seq[pos+1:]), splice_sites)
            for i, (changed, base) in enumerate(zip(changed_sites, base_sites)): per_site_changes[i].append(changed-base)
        for values, changes in zip(site_position_values, per_site_changes): values.append(sum(changes)/len(changes))
    weights = list(kernel) if kernel else [1.0]*(end-start)
    if len(weights) != end-start or any(w < 0 for w in weights) or sum(weights) <= 0: raise ValueError("footprint kernel must be nonnegative, nonzero, and match footprint")
    site_effects = tuple(sum(w*x for w, x in zip(weights, values))/sum(weights) for values in site_position_values)
    return {"untreated_propensity": base_total, "occlusion_importance": sum(site_effects), "site_untreated": base_sites, "site_occlusion": site_effects, "site_position_importance": tuple(tuple(x) for x in site_position_values), "contrast_kind": "model-derived-sequence-occlusion-not-native-treatment-output", "mapping": tuple((s.site_type, s.position, s.endpoint_direction) for s in splice_sites)}


# Sugimoto RNA/DNA nearest-neighbor delta-G37 parameters (kcal/mol).
RNA_DNA_DG37 = {"AA":-1.0,"AC":-2.1,"AG":-1.8,"AT":-0.9,"CA":-0.9,"CC":-2.1,"CG":-1.7,"CT":-0.9,"GA":-1.3,"GC":-2.7,"GG":-2.9,"GT":-1.1,"TA":-0.6,"TC":-1.5,"TG":-1.6,"TT":-0.2}
CHEMISTRY_PER_BASE = {"dna":0.0, "2ome":-0.20, "moe":-0.28, "lna":-0.45, "pmo":-0.10}


def hybridization_dg(sequence, oligo, chemistry="dna", monovalent_salt_molar=0.15):
    target, oligo = sequence.upper().replace("U", "T"), oligo.upper().replace("U", "T")
    if not target or not oligo: raise ValueError("target and oligo sequences are required")
    if set(target+oligo)-set("ACGT") or len(target) != len(oligo): raise ValueError("aligned canonical target/oligo of equal length required")
    if oligo != target.translate(str.maketrans("ACGT", "TGCA"))[::-1]: raise ValueError("oligo must be reverse-complementary to target footprint")
    key = chemistry.lower()
    if key not in CHEMISTRY_PER_BASE: raise ValueError("unregistered chemistry correction")
    if monovalent_salt_molar <= 0: raise ValueError("positive monovalent salt concentration required")
    return 2.2 + sum(RNA_DNA_DG37[target[i:i+2]] for i in range(len(target)-1)) + len(target)*CHEMISTRY_PER_BASE[key] + 0.114*max(0, len(target)-1)*math.log(monovalent_salt_molar)


# Complete ESEfinder 2.0 functional-SELEX matrices and canonical thresholds.
ESEFINDER={
 "SF2ASF":(((-1.14,.62,-1.58,1.32,-1.58,-1.58,.62),(1.37,-1.1,.73,.33,.94,-1.58,-1.58),(-.21,.17,.48,-1.58,.33,.99,-.11),(-1.58,-.5,-1.58,-1.13,-1.58,-1.13,.27)),1.956),
 "SC35":(((-.88,.09,-.06,-1.58,.09,-.41,-.06,.23),(-1.16,-1.58,.95,1.11,.56,.86,.32,-1.58),(.87,.45,-1.36,-1.58,-.33,-.05,-1.36,.68),(-1.18,-.2,.38,.88,-.2,-.86,.96,-1.58)),2.383),
 "SRp40":(((-.13,-1.58,1.28,-.33,.97,-.13,-1.58),(.56,.68,-1.12,1.24,-.77,.13,-.05),(-1.58,-.14,-1.33,-.48,-1.58,.44,.8),(.92,.37,.23,-1.14,.72,-1.58,-1.58)),2.670),
 "SRp55":(((-.66,.11,-.66,.11,-1.58,.61),(.39,-1.58,1.48,-1.58,-1.58,.98),(-1.58,.72,-1.58,.72,.21,-.79),(1.22,-1.58,-.07,-1.58,1.02,-1.58)),2.676)}
ESS_LIBRARY=(("hnRNPA1","TAGACA"),("hnRNPA1","TAGAGT"),("hnRNPA1","TAGGGW"),("hnRNPAB","ATAGCA"),("hnRNPH","GGCGG"),("hnRNPH","GGGTG"),("hnRNPH","TGTGGG"),("hnRNPH","TTGGGT"),("MBNL","YGCTKY"),("NOVA1","YCAY"),("PTBP1","CTCTCT"),("PTBP1","TCTT"),("CUGBP","TGCTG"),("YB1","CAACCACAA"),("FOX1_2","TGCATG"))
IUPAC={"A":"A","C":"C","G":"G","T":"T","U":"T","S":"GC","N":"ACGT","D":"AGT","Y":"CT","W":"AT","K":"GT"}


def motif_features(sequence):
    seq=sequence.upper().replace("U","T"); out={}; base_index={b:i for i,b in enumerate("ACGT")}
    for name,(matrix,threshold) in ESEFINDER.items():
        width=len(matrix[0]); scores=[sum(matrix[base_index[b]][j] for j,b in enumerate(seq[i:i+width])) for i in range(max(0,len(seq)-width+1))]
        out[f"motif_ESE_{name}_mean"]=sum(scores)/len(scores) if scores else 0.0
        out[f"motif_ESE_{name}_max"]=max(scores) if scores else 0.0
        out[f"motif_ESE_{name}_hit_density"]=sum(x>=threshold for x in scores)/len(scores) if scores else 0.0
    for name,pattern in ESS_LIBRARY:
        hits=sum(all(base in IUPAC[code] for base,code in zip(seq[i:i+len(pattern)],pattern)) for i in range(max(0,len(seq)-len(pattern)+1)))
        key=f"motif_ESS_{name}_density"; out[key]=out.get(key,0.0)-hits/max(1,len(seq)-len(pattern)+1)
    return out


def kmer_frequencies(sequence, k):
    seq = sequence.upper().replace("U", "T")
    if k < 1 or len(seq) < k: raise ValueError("k must fit sequence")
    total=len(seq)-k+1; found={}
    for i in range(total): found[seq[i:i+k]]=found.get(seq[i:i+k],0.0)+1.0/total
    return found


class FeatureBuilder:
    def __init__(self, cis_backends, accessibility_backend, k=(3,4,5)):
        self.cis_backends,self.accessibility_backend=tuple(cis_backends),accessibility_backend
        self.k_grid=(int(k),) if isinstance(k,int) else tuple(map(int,k))
    def build(self,item):
        target=item.sequence[item.footprint_start:item.footprint_end]; center=(item.footprint_start+item.footprint_end-1)/2
        features={"affinity_dg":hybridization_dg(target,item.oligo_sequence,item.chemistry,item.monovalent_salt_molar),"accessibility":self.accessibility_backend.accessibility(item.sequence,item.footprint_start,item.footprint_end),"length":float(item.footprint_end-item.footprint_start),"dose":item.dose,"chemistry_"+item.chemistry.lower():1.0}
        features.update(motif_features(target))
        for site in item.splice_sites: features[f"splice_distance_{site.site_type}_{'include' if site.endpoint_direction>0 else 'skip'}"]=site.endpoint_direction*(center-site.position)
        for k in self.k_grid:
            for mer,value in kmer_frequencies(target,k).items(): features[f"kmer_k{k}_"+mer]=value
        metadata={}
        for backend in self.cis_backends:
            occ=occlusion_importance(item.sequence,item.footprint_start,item.footprint_end,item.splice_sites,backend,item.footprint_kernel)
            for i,site in enumerate(item.splice_sites):
                stem=f"{backend.name}_{site.site_type}_{site.position}"; features[stem+"_untreated"]=occ["site_untreated"][i]; features[stem+"_occlusion"]=occ["site_occlusion"][i]
            features[f"{backend.name}_endpoint_occlusion"]=occ["occlusion_importance"]
            metadata[backend.name]={"contrast_kind":occ["contrast_kind"],"mapping":occ["mapping"]}
        return features,metadata


FEATURE_SETS={"affinity_dg":("affinity_dg",),"accessibility":("accessibility",),"motif_model":("motif_",),"splice_site_distance":("splice_distance_",),"kmer":("kmer_",),"cis_occlusion":("_untreated","_occlusion"),"nonstructural_stack":("affinity_dg","accessibility","motif_","splice_distance_","length","dose","chemistry_","kmer_","_untreated","_occlusion")}
