signals = {
    "indicators": [
        {
            "type": "sma",
            "name": "sma1",
            "args": 10
        },
        {
            "type": "sma",
            "name": "sma2",
            "args": 50
        },
        {
            "name": "sma3",
            "type": "sma",
            "args": 100
        },
        {
            "name": "rsi",
            "type": "rsi"
        }],
    "buy_signals": [
        {
            "name": "bsig1",
            "form": ("sma1", ">", "sma2")
        }
    ]
}
