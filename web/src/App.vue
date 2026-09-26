<template>
  <q-layout view="hHh Lpr lFf">
    <q-header class="site-header">
      <q-toolbar class="site-toolbar">
        <q-btn flat round dense icon="menu" class="menu-button" @click="drawerOpen = !drawerOpen" />
        <div class="brand">
          <div class="brand-mark">
            <span></span><span></span><span></span>
          </div>
          <div>
            <div class="brand-name">BioSTAR</div>
            <div class="brand-caption">BIOINFORMATICS RESOURCE</div>
          </div>
        </div>
        <q-space />
        <q-btn flat round icon="dark_mode" class="theme-button" @click="toggleDarkMode" />
      </q-toolbar>
    </q-header>

    <q-drawer v-model="drawerOpen" show-if-above :width="232" class="site-drawer">
      <div class="drawer-inner">
        <div class="drawer-title">TOOLS</div>

        <q-item clickable v-ripple to="/" exact class="nav-link">
          <q-item-section avatar><q-icon name="home" /></q-item-section>
          <q-item-section>Overview</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/sequences" class="nav-link">
          <q-item-section avatar><q-icon name="biotech" /></q-item-section>
          <q-item-section>Sequences</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/proteins" class="nav-link">
          <q-item-section avatar><q-icon name="science" /></q-item-section>
          <q-item-section>Proteins</q-item-section>
        </q-item>

        <q-item clickable v-ripple to="/mutations" class="nav-link">
          <q-item-section avatar><q-icon name="compare_arrows" /></q-item-section>
          <q-item-section>Mutations</q-item-section>
        </q-item>

        <div class="drawer-divider"></div>
        <div class="drawer-title">PLATFORM</div>

        <q-item clickable v-ripple to="/settings" class="nav-link">
          <q-item-section avatar><q-icon name="settings" /></q-item-section>
          <q-item-section>Settings</q-item-section>
        </q-item>

        <q-space />
        <div class="drawer-footer">
          <div>BioSTAR API</div>
          <span>Scientific computing tools</span>
        </div>
      </div>
    </q-drawer>

    <q-page-container>
      <router-view />
    </q-page-container>
  </q-layout>
</template>

<script setup lang="ts">
import { ref } from "vue";
import { Dark } from "quasar";

const drawerOpen = ref<boolean>(true);

const toggleDarkMode = (): void => {
  Dark.toggle();
};
</script>

<style scoped>
.site-header {
  background: #123d59;
  color: white;
}

.site-toolbar {
  min-height: 68px;
  padding: 0 22px;
}

.menu-button,
.theme-button {
  color: rgb(255 255 255 / 78%);
}

.brand {
  display: flex;
  align-items: center;
  gap: 11px;
  margin-left: 12px;
}

.brand-mark {
  display: flex;
  align-items: flex-end;
  gap: 2px;
  width: 30px;
  height: 30px;
}

.brand-mark span {
  width: 7px;
  border-radius: 5px 5px 2px 2px;
  background: #5bc6c0;
  transform: skewY(-22deg);
}

.brand-mark span:nth-child(1) { height: 15px; }
.brand-mark span:nth-child(2) { height: 23px; }
.brand-mark span:nth-child(3) { height: 30px; }

.brand-name {
  font-size: 1.15rem;
  font-weight: 800;
  letter-spacing: 0.06em;
}

.brand-caption {
  margin-top: 1px;
  color: rgb(255 255 255 / 58%);
  font-size: 0.56rem;
  letter-spacing: 0.13em;
}

.site-drawer {
  background: #fff;
  border-right: 1px solid #e4eaee;
}

.drawer-inner {
  display: flex;
  min-height: 100%;
  flex-direction: column;
  padding: 28px 12px 18px;
}

.drawer-title {
  padding: 0 13px 8px;
  color: #82929c;
  font-size: 0.65rem;
  font-weight: 800;
  letter-spacing: 0.14em;
}

.nav-link {
  min-height: 44px;
  margin: 2px 0;
  border-radius: 5px;
  color: #435963;
  font-size: 0.9rem;
}

.nav-link .q-icon {
  color: #82949d;
}

.nav-link.q-router-link--active {
  background: #e8f4f5;
  color: #126d79;
  font-weight: 700;
}

.nav-link.q-router-link--active .q-icon {
  color: #168391;
}

.drawer-divider {
  height: 1px;
  margin: 18px 12px;
  background: #e4eaee;
}

.drawer-footer {
  margin: 18px 12px 0;
  padding-top: 14px;
  border-top: 1px solid #e4eaee;
  color: #536872;
  font-size: 0.75rem;
  font-weight: 700;
}

.drawer-footer span {
  display: block;
  margin-top: 3px;
  color: #91a0a8;
  font-size: 0.68rem;
  font-weight: 400;
}
</style>
