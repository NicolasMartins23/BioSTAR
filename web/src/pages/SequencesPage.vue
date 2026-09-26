<template>
  <q-page padding>
    <div class="page-shell">
      <div class="text-overline text-primary">Analysis</div>
      <div class="text-h4 text-weight-bold">Sequences</div>
      <div class="text-body1 text-grey-6 q-mt-sm">
        Convert DNA and RNA sequences using the BioSTAR API.
      </div>

      <q-card flat bordered class="q-mt-xl">
        <q-card-section>
          <div class="text-h6">Sequence conversion</div>
          <div class="text-body2 text-grey-6 q-mt-xs">
            Enter a single DNA or RNA sequence. FASTA input is also supported.
          </div>

          <q-select
            v-model="conversion"
            :options="conversionOptions"
            emit-value
            map-options
            outlined
            label="Conversion"
            class="q-mt-lg"
          />

          <q-input
            v-model="sequence"
            outlined
            type="textarea"
            autogrow
            label="Input sequence"
            hint="Whitespace is ignored. Maximum length: 1,000 characters."
            class="q-mt-md"
            :disable="isLoading"
            @keydown.ctrl.enter="convert"
          />

          <div class="row justify-end q-mt-md">
            <q-btn
              unelevated
              color="primary"
              icon="transform"
              label="Convert"
              :loading="isLoading"
              :disable="!sequence.trim()"
              @click="convert"
            />
          </div>
        </q-card-section>
      </q-card>

      <q-card
        v-if="errorMessage"
        flat
        bordered
        class="q-mt-md"
      >
        <q-card-section>
          <q-banner rounded class="bg-negative text-white">
            {{ errorMessage }}
          </q-banner>
        </q-card-section>
      </q-card>

      <q-card
        v-if="result"
        flat
        bordered
        class="q-mt-md"
      >
        <q-card-section>
          <div class="row items-center justify-between">
            <div>
              <div class="text-h6">Result</div>
              <div class="text-caption text-grey-6">{{ conversionLabel }}</div>
            </div>
            <q-btn
              flat
              round
              icon="content_copy"
              aria-label="Copy result"
              @click="copyResult"
            >
              <q-tooltip>Copy result</q-tooltip>
            </q-btn>
          </div>

          <q-input
            :model-value="result"
            readonly
            outlined
            type="textarea"
            autogrow
            class="q-mt-md"
          />
        </q-card-section>
      </q-card>
    </div>
  </q-page>
</template>

<script setup lang="ts">
import { computed, ref } from "vue";
import { Notify } from "quasar";
import type { APIResponse, SequenceResponse } from "../services/api";

import {
  convertDnaToProtein,
  convertDnaToRna,
  convertRnaToDna,
  convertRnaToProtein,
} from "../services/api";

type Conversion = "dna-rna" | "dna-protein" | "rna-protein" | "rna-dna";

const conversion = ref<Conversion>("dna-rna");
const sequence = ref<string>("");
const result = ref<string>("");
const errorMessage = ref<string>("");
const isLoading = ref<boolean>(false);

const conversionOptions: Array<{ label: string; value: Conversion }> = [
  { label: "DNA → RNA", value: "dna-rna" },
  { label: "DNA → Protein", value: "dna-protein" },
  { label: "RNA → Protein", value: "rna-protein" },
  { label: "RNA → DNA", value: "rna-dna" },
];

const conversionLabel = computed<string>(() => {
  const option = conversionOptions.find(
    (item) => item.value === conversion.value,
  );

  return option?.label ?? conversion.value;
});

const convert = async (): Promise<void> => {
  const input: string = sequence.value.trim();

  if (!input) {
    return;
  }

  isLoading.value = true;
  errorMessage.value = "";
  result.value = "";

  try {
    let response: APIResponse<SequenceResponse>;

    if (conversion.value === "dna-rna") {
      response = await convertDnaToRna(input);
    } else if (conversion.value === "dna-protein") {
      response = await convertDnaToProtein(input);
    } else if (conversion.value === "rna-protein") {
      response = await convertRnaToProtein(input);
    } else {
      response = await convertRnaToDna(input);
    }

    if (response.data === null) {
      errorMessage.value = response.message?.message ?? "Sequence conversion returned no data.";
      return;
    }

    result.value = response.data.sequence;
  } catch (error: unknown) {
    errorMessage.value = error instanceof Error
      ? error.message
      : "Sequence conversion failed.";
  } finally {
    isLoading.value = false;
  }
};

const copyResult = async (): Promise<void> => {
  if (!result.value) {
    return;
  }

  try {
    await navigator.clipboard.writeText(result.value);
    Notify.create({
      type: "positive",
      message: "Result copied to clipboard.",
    });
  } catch {
    Notify.create({
      type: "negative",
      message: "Unable to copy the result.",
    });
  }
};
</script>

<style scoped>
.page-shell {
  max-width: 1200px;
  margin: 0 auto;
}
</style>
