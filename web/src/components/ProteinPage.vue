<template>
<div>
  <div class="q-mb-lg"><div class="text-overline text-primary bio-mono">PROTEIN</div><div class="text-h4 brand-title q-mt-sm">Protein Analysis</div><div class="text-body2 text-grey-7 q-mt-sm">Select the calculations you need, or run the complete analysis.</div></div>
  <q-card class="bio-card q-pa-lg">
    <q-input v-model="sequence" outlined type="textarea" autogrow :rows="9" label="Protein sequence" hint="Protein sequence or a single FASTA record" spellcheck="false" class="sequence-input" />
    <div class="row items-center q-mt-lg q-gutter-md"><q-checkbox v-model="fullResults" label="Full test results" /><q-input v-model.number="chargePH" type="number" outlined dense label="Charge at pH" style="width:160px" :disable="fullResults" /></div>
    <div class="row q-col-gutter-sm q-mt-md"><div v-for="option in options" :key="option.key" class="col-12 col-sm-6 col-md-4"><q-checkbox v-model="selected[option.key]" :label="option.label" :disable="fullResults" /></div></div>
    <q-btn color="primary" no-caps class="q-mt-lg full-width" label="Analyze protein" :loading="loading" :disable="!sequence.trim()" @click="run" />
  </q-card>
  <q-card class="bio-card q-pa-lg q-mt-lg"><div class="text-overline text-grey-6 bio-mono">RESULT</div>
    <q-banner v-if="error" rounded class="bg-red-1 text-negative q-mt-md">{{ error }}</q-banner>
    <ResultView v-else-if="result" :result="result" class="q-mt-md" />
    <div v-else class="text-grey-6 q-py-xl text-center">Run an analysis to see the result.</div>
  </q-card>
</div>
</template>
<script setup lang="ts">
import { reactive, ref } from "vue"; import { Notify } from "quasar"; import ResultView from "./ResultView.vue"; import { apiPost } from "../services/api";
type AnalysisKey="aminoacids_count"|"isoelectric_point"|"aromaticity"|"secondary_structure_propensity"|"molecular_weight"|"hydrophobic_index"|"composition_ratio"|"extinction_coefficient";
const sequence=ref("");const fullResults=ref(false);const chargePH=ref<number|null>(7);const loading=ref(false);const error=ref("");const result=ref<Record<string,unknown>|null>(null);
const options:Array<{key:AnalysisKey;label:string}>=[{key:"aminoacids_count",label:"Amino-acid count"},{key:"isoelectric_point",label:"Isoelectric point"},{key:"aromaticity",label:"Aromaticity"},{key:"secondary_structure_propensity",label:"Secondary structure"},{key:"molecular_weight",label:"Molecular weight"},{key:"hydrophobic_index",label:"Hydrophobic index"},{key:"composition_ratio",label:"Composition ratio"},{key:"extinction_coefficient",label:"Extinction coefficient"}];
const selected=reactive<Record<AnalysisKey,boolean>>({aminoacids_count:false,isoelectric_point:false,aromaticity:false,secondary_structure_propensity:false,molecular_weight:false,hydrophobic_index:false,composition_ratio:false,extinction_coefficient:false});
async function run():Promise<void>{loading.value=true;error.value="";result.value=null;const payload:Record<string,unknown>={sequence:sequence.value,get_full_test_results:fullResults.value,get_aminoacids_count:selected.aminoacids_count,get_isoelectric_point:selected.isoelectric_point,get_aromaticity:selected.aromaticity,get_secondary_structure_propensity:selected.secondary_structure_propensity,get_molecular_weight:selected.molecular_weight,get_hydrophobic_index:selected.hydrophobic_index,get_composition_ratio:selected.composition_ratio,get_extinction_coefficient:selected.extinction_coefficient};if(!fullResults.value&&chargePH.value!==null)payload.get_charge_at_pH=chargePH.value;try{result.value=await apiPost("/api/protein",payload);Notify.create({type:"positive",message:"Protein analysis completed"});}catch(cause){error.value=cause instanceof Error?cause.message:"Request failed";}finally{loading.value=false;}}
</script>
