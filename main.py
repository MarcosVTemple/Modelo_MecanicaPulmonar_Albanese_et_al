from modelos.main_mec_pulm import MecanicaPulmonar

if __name__ == "__main__":
    mec_pulm_obj = MecanicaPulmonar()
    mec_pulm_obj.run_mecanica_pulmonar()
    mec_pulm_obj.plot_mecanica_pulmonar()
