using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class InfoBoxDisplay : MonoBehaviour
{
    Text text;
    private void Awake()
    {
        text = GetComponent<Text>();
    }
    void Update()
    {
        text.text = "B_mean= " + SimulationManager.Bmean.ToString("F3") + "\n";
        text.text += "N_spots= " + SimulationManager.Nspots.ToString() + "\n";
        text.text += "MaxConcentration= " + SimulationManager.MaxConcentration.ToString() + "\n";
    }
}

