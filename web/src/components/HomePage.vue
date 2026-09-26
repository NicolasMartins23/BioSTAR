<template>
<div>
  <div class="row q-col-gutter-lg">
    <div class="col-12 col-lg-7"><q-card class="bio-card q-pa-xl">
      <div class="text-overline text-primary bio-mono">TARGETED BIOINFORMATICS</div>
      <h1 class="brand-title text-h2 q-my-md">BioSTAR</h1>
      <p class="text-body1 text-grey-7">Bioinformatics software for targeted DNA, RNA and protein analysis. Run sequence conversions and focused protein analysis through the BioSTAR API.</p>
      <div class="q-gutter-sm q-mt-lg"><q-btn color="primary" no-caps label="DNA → Protein" @click="$emit('select','dna-protein')" /><q-btn outline color="primary" no-caps label="Protein Analysis" @click="$emit('select','protein')" /></div>
    </q-card></div>
    <div class="col-12 col-lg-5"><q-card class="bio-card q-pa-lg full-height">
      <div class="text-overline text-grey-6 bio-mono">QUICK START</div><div class="text-h6 q-mt-sm">Try a conversion</div>
      <div class="text-body2 text-grey-7 q-mt-sm">Short interactive conversions use the API GET endpoints and accept DNA/RNA sequences or a single FASTA record.</div>
      <div class="q-gutter-sm q-mt-lg"><q-btn flat color="primary" no-caps label="DNA → RNA" @click="$emit('select','dna-rna')" /><q-btn flat color="primary" no-caps label="RNA → Protein" @click="$emit('select','rna-protein')" /></div>
    </q-card></div>
  </div>
  <div class="row q-col-gutter-lg q-mt-lg">
    <div v-for="item in tools" :key="item.title" class="col-12 col-md-4"><q-card class="bio-card q-pa-lg full-height">
      <q-icon :name="item.icon" size="28px" color="primary" /><div class="text-h6 q-mt-md">{{ item.title }}</div><div class="text-body2 text-grey-7 q-mt-sm">{{ item.description }}</div>
      <q-btn flat no-caps color="primary" class="q-px-none q-mt-md" :label="item.action" @click="$emit('select',item.tool)" />
    </q-card></div>
  </div>
</div>
</template>
<script setup lang="ts">
defineEmits<{select:[tool:"dna-rna"|"dna-protein"|"rna-protein"|"rna-dna"|"protein"|"mutation"]}>();
const tools=[
 {icon:"swap_horiz",title:"Sequence Conversion",description:"Convert DNA and RNA sequences or translate DNA/RNA into protein.",action:"Open conversions",tool:"dna-rna" as const},
 {icon:"biotech",title:"Protein Analysis",description:"Request targeted protein metrics or a full analysis result.",action:"Analyze protein",tool:"protein" as const},
 {icon:"compare_arrows",title:"Mutation Compare",description:"Compare two coding DNA sequences and inspect detected mutations.",action:"Compare sequences",tool:"mutation" as const},
];
</script>
