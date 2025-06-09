#!/usr/bin/env python3

import math
import json
import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from pathlib import Path

def write_cdc_custom_file(filename, fields, fov):
    fov_deg = max(fov) * 60.0
    with open(filename, "w") as f:
        for ra, dec, name in fields:
            dec_fmt = '+' + dec if not dec.strip().startswith('-') else dec
            f.write(f"{name},{ra},{dec_fmt},1,16,0.00,{fov_deg},30975,{name}\n")

def create_target(index, ra, dec, name, plan_name, exposure, count):
    target = ET.Element(f"Target{index}")
    target.attrib = {
        "PA": "-",
        "RA": ra,
        "Dec": dec,
        "Path": "",
        "Plan": plan_name,
        "Skip": "False",
        "Delay": "1",
        "EndSet": "True",
        "EndTime": "",
        "Preview": "False",
        "FlatBinX": "0",
        "FlatBinY": "0",
        "FlatGain": "0",
        "FullOnly": "False",
        "DarkNight": "False",
        "FlatCount": "0",
        "FlatFstop": "",
        "Magnitude": "-",
        "StartRise": "True",
        "StartTime": "",
        "FlatOffset": "0",
        "ObjectName": name,
        "RepeatDone": "0",
        "ScriptArgs": "",
        "EndMeridian": "-9999",
        "FlatFilters": "",
        "HFM_Enabled": "False",
        "RepeatCount": "1",
        "UpdateCoord": "False",
        "AutofocusTemp": "False",
        "MoonAvoidance": "True",
        "SolarTracking": "False",
        "StartMeridian": "-9999",
        "PreviewExposure": "1",
        "InplaceAutofocus": "False",
        "AstrometryPointing": "True",
        "MandatoryStartTime": "False",
        "NoAutoguidingChange": "False"
    }

    init_script = ET.SubElement(target, "InitScript")
    init_script.attrib = {
        "InitScript": "False",
        "InitScriptArgs": "",
        "InitScriptName": "",
        "InitScriptPath": "/home/tstibor/.config/ccdciel/"
    }

    plan = ET.SubElement(target, "Plan")
    plan.attrib = {"Name": plan_name, "StepNum": "1"}

    steps = ET.SubElement(plan, "Steps")

    step1 = ET.SubElement(steps, "Step1")
    step1.attrib = {
        "Done": "0",
        "Gain": "1",
        "Type": "0",
        "Count": count,
        "Fstop": "",
        "Dither": "True",
        "Filter": "No change",
        "Offset": "0",
        "Binning": "1x1",
        "Exposure": exposure,
        "Autofocus": "False",
        "FrameType": "Light",
        "ScriptArgs": "",
        "ScriptName": "",
        "ScriptPath": "",
        "StackCount": "1",
        "SwitchName": "",
        "Description": f"{count}x{exposure}sec",
        "DitherCount": "1",
        "RefExposure": "",
        "SwitchValue": "",
        "AutofocusCount": "0",
        "AutofocusStart": "False",
        "SwitchNickname": ""
    }

    return target

def repeat_fields_in_order(group_fields, repeat_count):
    return [field for _ in range(repeat_count) for field in group_fields]

def sexagesimal_to_decimal_ra(ra_str):
    h, m, s = ra_str.replace("h", " ").replace("m", " ").replace("s", "").split()
    return round((float(h) + float(m) / 60 + float(s) / 3600) * 15, 6)

def sexagesimal_to_decimal_dec(dec_str):
    sign = -1 if dec_str.strip().startswith("-") else 1
    d = dec_str.replace("-", "").replace("+", "")
    d, m, s = d.replace("d", " ").replace("m", " ").replace("s", "").split()
    return round(sign * (float(d) + float(m) / 60 + float(s) / 3600), 6)

def decimal_to_sexagesimal_ra(ra_deg):
    total_seconds = ra_deg / 15 * 3600
    h = int(total_seconds // 3600)
    m = int((total_seconds % 3600) // 60)
    s = total_seconds % 60
    return f"{h:02d}h{m:02d}m{s:05.2f}s"

def generate_ra_dec(base_ra_sex, base_dec_sex, text, fov_width_deg, num_tiles):
    base_ra_deg = sexagesimal_to_decimal_ra(base_ra_sex)
    base_dec_deg = sexagesimal_to_decimal_dec(base_dec_sex)
    ra_step = fov_width_deg / math.cos(math.radians(base_dec_deg))

    ra_list = []
    for i in range(num_tiles):  # Number of horizontal tiles.
        new_ra_deg = base_ra_deg + i * ra_step
        new_ra_sex = decimal_to_sexagesimal_ra(new_ra_deg)
        ra_list.append((new_ra_sex, base_dec_sex, f"{text} Tile {i+1}"))

    return ra_list

def calc_fov_deg(focal_len_mm, sensor_pixels_x, sensor_pixels_y, sensor_um_x, sensor_um_y):
    pixel_scale_x = (206.265 * sensor_um_x) / focal_len_mm
    pixel_scale_y = (206.265 * sensor_um_y) / focal_len_mm

    fov_x_deg = sensor_pixels_x * pixel_scale_x / 3600.0
    fov_y_deg = sensor_pixels_y * pixel_scale_y / 3600.0

    return (fov_x_deg, fov_y_deg)

def write_obslist_file(filename, fields, list_title="Generated Observing List"):
    with open(filename, "w") as f:
        f.write(f"{list_title}\n")
        for ra, dec, name in fields:
            ra_deg = f"{sexagesimal_to_decimal_ra(ra):10.6f}"
            dec_deg = f"{sexagesimal_to_decimal_dec(dec):10.6f}"
            f.write(
                f"{name:<32}"       # object name
                f"{ra_deg:>10}"     # RA decimal
                f"{dec_deg:>10}"    # Dec decimal
                f"{name:<32}"       # object name again for label/update
                f"Field generated from {name}\n"  # optional comment
            )

def generate_all_outputs(all_fields, base_name, fov, plan_name, args):

    config = ET.Element("CONFIG", Version="5", ListName="Fields", TargetNum=str(len(all_fields)), RepeatCount="1")
    targets = ET.SubElement(config, "Targets", RepeatDone="0", ResetRepeat="True", IgnoreRestart="False")

    for idx, (ra, dec, name) in enumerate(all_fields, 1):
        targets.append(create_target(idx, ra, dec, name, plan_name, args.exposure, args.count))

    ET.SubElement(config, "Startup", {
        "Unpark": "True",
        "SeqStop": "True",
        "SeqStart": "True",
        "RunScript": "False",
        "CoolCamera": "True",
        "SeqStopTwilight": "True",
        "SeqStartTwilight": "True",
        "StartScript": ""
    })

    ET.SubElement(config, "Termination", {
        "Park": "True",
        "CloseDome": "False",
        "EndScript": "",
        "RunScript": "False",
        "WarmCamera": "True",
        "ErrorScript": "",
        "StopTracking": "True",
        "ErrorRunScript": "False"
    })

    # Write ccdciel targets sequence file.
    xml_str = ET.tostring(config, encoding="utf-8")
    pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="  ")
    Path(f"{base_name}.targets").write_text(pretty_xml)
    print(f"✅ Created {base_name}.targets")

    # Write skychart observation list.
    write_obslist_file(f"{base_name}.txt", all_fields)
    print(f"✅ Created {base_name}.txt")

    # Write custom object skychart file.
    write_cdc_custom_file(f"{base_name}.custom", all_fields, fov)
    print(f"✅ Created {base_name}.custom")

def load_constellations_fields(json_path):
    with open(json_path, "r") as f:
        return json.load(f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate observation plans from constellation fields.")
    parser.add_argument(
        "json_file",
        nargs="?",
        default="fields.json",
        help="Path to the constellation fields JSON file (default: fields.json)"
    )
    parser.add_argument(
        "--focal-length", type=float, default=600.0,
        help="Telescope focal length in mm (default: 600.0)"
    )
    parser.add_argument(
        "--sensor-width-px", type=int, default=1391,
        help="Sensor width in pixels (default: 1391)"
    )
    parser.add_argument(
        "--sensor-height-px", type=int, default=1039,
        help="Sensor height in pixels (default: 1039)"
    )
    parser.add_argument(
        "--pixel-size-um", type=float, default=6.45,
        help="Sensor pixel size in microns (assumes square pixels, default: 6.45)"
    )
    parser.add_argument(
        "--repeat", type=int, default=3,
        help="Repetition of each field (default: 3)"
    )
    parser.add_argument(
        "--exposure", type=str, default="60",
        help="Exposure time in seconds (default: 60)"
    )
    parser.add_argument(
        "--count", type=str, default="5",
        help="Number of exposures (default: 5)"
    )
    
    args = parser.parse_args()
    plan_name = f"{args.count}x{args.exposure}sec"
    constellations_fields = load_constellations_fields(args.json_file)

    fov = calc_fov_deg(
        args.focal_length,
        args.sensor_width_px,
        args.sensor_height_px,
        args.pixel_size_um,
        args.pixel_size_um
    )
    fov_x_with_extend = fov[0] # + 0.1

    all_fields = []
    for constellation, fields in constellations_fields.items():
        for field in fields:
            tiles = generate_ra_dec(field[0], field[1], field[2], fov_x_with_extend, 3)
            all_fields.extend(tiles * args.repeat)

    generate_all_outputs(all_fields, "fields", fov, plan_name, args)
