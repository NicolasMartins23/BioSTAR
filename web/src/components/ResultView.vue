<template><div class="row q-col-gutter-sm"><div v-for="[key,value] in entries" :key="key" class="col-12 col-md-6"><q-card flat bordered class="q-pa-md"><div class="text-caption text-grey-6 text-uppercase bio-mono">{{formatKey(key)}}</div><pre v-if="typeof value==='object'" class="q-mt-sm result-pre">{{JSON.stringify(value,null,2)}}</pre><div v-else class="text-subtitle2 text-primary q-mt-sm">{{String(value)}}</div></q-card></div></div></template>
<script setup lang="ts">
import { computed } from "vue"; const props=defineProps<{result:Record<string,unknown>}>(); const entries=computed(()=>Object.entries(props.result)); function formatKey(key:string):string{return key.replaceAll("_"," ");}
</script>
