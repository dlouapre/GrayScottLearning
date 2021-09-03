using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class TimeDisplay : MonoBehaviour
{
    Text text;
    private void Awake()
    {
        text = GetComponent<Text>();
    }
    void Update()
    {
        text.text = "t=" + SimulationManager.time.ToString("F0");
    }
}

