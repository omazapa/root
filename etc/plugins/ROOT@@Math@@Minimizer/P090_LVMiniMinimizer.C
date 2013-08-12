void P090_LVMiniMinimizer()
{
   gPluginMgr->AddHandler("ROOT::Math::Minimizer", "LVMini", "ROOT::LVMini::LVMiniMinimizer",
      "LVMini", "LVMiniMinimizer(const char *)");
}
