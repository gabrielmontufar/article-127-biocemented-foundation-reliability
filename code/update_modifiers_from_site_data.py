from modifier_site_update_synthetic_example import (
    calibrate_kh_from_plate_load,
    calibrate_kz_from_cpt_profile,
    modifier_update_and_sensitivity,
    write_readme,
)


if __name__ == "__main__":
    kh = calibrate_kh_from_plate_load()
    kz = calibrate_kz_from_cpt_profile()
    base = modifier_update_and_sensitivity(kh["kH_fit"], kz["kz_fit"])
    write_readme()
    print({"kH": kh, "kz": kz, "baseline": base})
