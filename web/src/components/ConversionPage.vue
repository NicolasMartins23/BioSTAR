<template>
<div>
  <div class="q-mb-lg"><div class="text-overline text-primary bio-mono">SEQUENCE CONVERSION</div><div class="text-h4 brand-title q-mt-sm">{{ title }}</div><div class="text-body2 text-grey-7 q-mt-sm">Maximum interactive sequence length: 1,000 nucleotides.</div></div>
  <q-card class="bio-card q-pa-lg">
    <q-input v-model="sequence" outlined type="textarea" autogrow :rows="9" label="Sequence" hint="DNA / RNA / single FASTA record" spellcheck="false" class="sequence-input" />
    <div class="row justify-between items-center q-mt-md"><div class="text-caption text-grey-6 bio-mono">{{ normalizedLength.toLocaleString() }} characters</div><q-btn color="primary" no-caps label="Run conversion" :loading="loading" :disable="!sequence.trim()" @click="run" /></div>
  </q-card>
  <q-card class="bio-card q-pa-lg q-mt-lg"><div class="text-overline text-grey-6 bio-mono">RESULT</div>
    <q-banner v-if="error" rounded class="bg-red-1 text-negative q-mt-md">{{ error }}</q-banner>
    <div v-else-if="result" class="q-mt-md"><q-input :model-value="String(result.sequence ?? '')" readonly outlined type="textarea" autogrow class="sequence-input" /></div>
    <div v-else class="text-grey-6 q-py-xl text-center">Run a conversion to see the result.</div>
  </q-card>
</div>
</template>
<script setup lang="ts">
import { computed, ref } from "vue";
import { Notify } from "quasar";
import { apiGet } from "../services/api";
type ConversionTool="dna-rna"|"dna-protein"|"rna-protein"|"rna-dna";
const props=defineProps<{tool:ConversionTool}>();
const endpointMap:Record<ConversionTool,string>={"dna-rna":"/api/dna-rna","dna-protein":"/api/dna-protein","rna-protein":"/api/rna-protein","rna-dna":"/api/rna-dna"};
const titleMap:Record<ConversionTool,string>={"dna-rna":"DNA → RNA","dna-protein":"DNA → Protein","rna-protein":"RNA → Protein","rna-dna":"RNA → DNA"};
const sequence=ref(""); const result=ref<Record<string,unknown>|null>(null); const error=ref(""); const loading=ref(false);
const title=computed(()=>titleMap[props.tool]); const normalizedLength=computed(()=>sequence.value.replace(/^>[^\n]*\n?/,"").replace(/\s/g,"").length);
async function run():Promise<void>{loading.value=true;error.value="";result.value=null;try{result.value=await apiGet(endpointMap[props.tool],{sequence:sequence.value});Notify.create({type:"positive",message:"Conversion completed"});}catch(cause){error.value=cause instanceof Error?cause.message:"Request failed";}finally{loading.value=false;}}
</script>
